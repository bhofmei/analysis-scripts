import sys, math, glob, multiprocessing, subprocess, os
from bioFiles import *

# Usage: python region_methylation_pe.py [-r] [-k] [-p=num_proc] [-m=meth_types] [-o=out_id] <region_file> <allc_path> <sample_name> [sample_name]*

NUMPROC=1

def processInputs( dmrFileStr, allCPath, lineNamesAr, readCounts, isPickle, outID, numProc, mTypes ):
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	print( 'Methylation types: {:s}\nSamples included: {:s}\nNumber of DMRs: {:d}\nRead counts: {:s}\nUse pickle: {:s}\n'.format( '  '.join(mTypes), '  '.join(lineNamesAr),dmrCount, str(readCounts), str(isPickle) ) )
	info = '#from_script: region_methylation_pe.py; meth_types: {:s}; read_counts: {:s}; region_file: {:s}'.format( ','.join(mTypes), str(readCounts), os.path.basename( dmrFileStr ) )
	
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processSample, args=(allCPath, dmrDict, sampleName, mTypes, readCounts, isPickle) ) for sampleName in lineNamesAr ]
	suc = [ p.get() for p in results ]
	
	outFileStr = outID + '_region_meth.tsv'
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, suc, info, readCounts )
	
	print( 'Done.' )


def readDMRFile( dmrFileStr ):
	'''
		return dictionary of subset (DMR) regions
		{chr1:[(start,end),(start2,end2)],chr2:[(start,end)],...}
	'''
	dmrFile = open( dmrFileStr, 'r' )
	dmrDict = {}
	dmrCount = 0
	
	for line in dmrFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split()
		chrm = lineAr[0]
		start = int( lineAr[1] )
		end = int( lineAr[2] )
		if dmrDict.get(chrm) == None:
			dmrDict[chrm] = []
		dmrDict[chrm] += [(start, end)]
		dmrCount += 1
	dmrFile.close()
	return dmrDict, dmrCount


def processSample( allCPath, dmrDict, sampleName, mTypes, readCounts, pickle ):
	# get allc dict
	allcFileAr = [ 'allc_'+sampleName+'.tsv', 'allc_'+sampleName+'.tsv.gz', sampleName+'.tsv', sampleName+'.tsv.gz']
	allcFileStr = None
	for a in allcFileAr:
		t = os.path.normpath( allCPath + '/' + a )
		if os.path.exists( t ):
			allcFileStr = t
			break
	if allcFileStr == None:
		print( 'WARNING: allc file for', sampleName, 'not found...ignoring' )
		return ''
	# get allc dict
	allcFile = FileAllC_full( allcFileStr )
	allcDict = allcFile.getAllCDict(isPickle=pickle)
	
	# loop through mTypes
	outStr = ''
	for m in mTypes:
		outStr += processSampleMethType( allcDict[m], dmrDict, readCounts, '{:s}\t{:s}'.format( m, sampleName) )
	# end for m
	del allcDict
	return outStr
	
def processSampleMethType( allcDict, dmrDict, readCounts, info ):
	outStr = ''
	# loop through chrms
	for chrm in dmrDict.keys():
		chrmDict = allcDict.get(chrm)
		if chrmDict == None:
			print( 'WARNING: chrm', chrm, 'not found in allC file' )
			return ''
		outStr += processSampleChrm( chrmDict, dmrDict[ chrm ], readCounts, info + '\t{:s}'.format( chrm ) )
	return outStr
		
def processSampleChrm( allcDict, regionAr, readCounts, info ):
	outStr = ''
	# loop through regions
	for start, end in regionAr:
		calc = calculateRegionC( start, end, allcDict, readCounts )
		outStr += info + '\t' + calc + '\n'
	return outStr

def calculateRegionC( start, end, allcDict, readCounts ):
	
	mC = 0
	tC = 0
	nC = 0
	for pos in range(start, end+1):
		tup = allcDict.get(pos)
		if tup != None:
			mC += tup[0]
			tC += tup[1]
			nC += 1
	# end for pos
	wMeth = (0.0 if tC == 0 else float(mC) / float(tC))
	outStr = '{:d}\t{:d}\t{:.6f}'.format( start, end, wMeth )
	if readCounts:
		outStr += '\t{:d}\t{:d}\t{:d}'.format( mC, tC, nC )
	return outStr

def writeOutput( outFileStr, outMat, info, readCounts ):
	headerAr = [ 'meth_type', 'sample', 'chrm', 'start', 'end', 'wei_meth' ]
	if readCounts:
		headerAr += ['meth_reads','total_reads','num_c']
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' + '\t'.join( headerAr ) + '\n' )
	for text in outMat:
		outFile.write( text )
	outFile.close()

def parseInputs( argv ):
	
	readCounts = False
	isPickle = True
	numProc = NUMPROC
	mTypes = ['C', 'CG', 'CHG', 'CHH']
	outID = 'regions_out'
	startInd = 0
	
	for i in range(min(5, len(argv)-3)):
		if argv[i] == '-r':
			readCounts = True
			startInd += 1
		elif argv[i] == '-k':
			isPickle = False
			startInd += 1
		elif argv[i].startswith('-m'):
			mTypes = argv[i][3:].upper().split(',')
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	
	dmrFileStr = argv[startInd]
	allCPath = argv[startInd+1]
	if os.path.isdir( allCPath ) == False:
		print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allCPath ) )
		exit()
		
	lineNamesAr = []
	for i in range(startInd+2, len(argv)):
		lineNamesAr += [ argv[i] ]
	processInputs( dmrFileStr, allCPath, lineNamesAr, readCounts, isPickle, outID, numProc, mTypes )
		

def printHelp():
	print( 'Usage: python region_methylation_pe.py [-r] [-k] [-p=num_proc] [-m=meth_types] [-o=out_id] <region_file> <allc_path> <sample_name> [sample_name]*' )
	

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
