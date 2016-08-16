import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: python dmr_gen_counts_pe.py [-o=outID] [-m=methTypes] [-p=numProc] [-v=minCov] <dmrFile> <allcPath>  <sample1> <sample2> [sampleN]*
# does the between-generation computations of methylation and outputs to a file
# that will be used by R; input order of samples matters
# allc files should have all chrms together
COV=None

def processInputs( dmrFileStr, allcPath, sampleNamesAr, outID, methType, numProc, minCov ):
	print( 'DMR file:', os.path.basename( dmrFileStr ) )
	print( 'AllC path:', allcPath )
	print( 'Samples:', ', '.join(sampleNamesAr) )
	print( 'Methylation type:',  methType )
	print( 'Minimum coverage:', minCov )
	info = '#from_script: dmr_gen_counts_pe.py; dmr_file: {:s}; samples: {:s}; methyl_type: {:s}; min_cov: {:d}'.format( os.path.basename( dmrFileStr ), ','.join( sampleNamesAr ), methType, minCov )
	
	# read DMR file
	print( 'Reading DMR file' )
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	print( ' analyzing', dmrCount, 'DMRs' )
	
	# loop through adjacent generations
	outMat = []
	mDict1 = None
	mDict2 = {}
	pool = multiprocessing.Pool( processes=numProc )
	for i in range(len(sampleNamesAr)-1):
		print( 'Comparing', sampleNamesAr[i], 'and', sampleNamesAr[i+1] )
		comparisonName = '{:s}-{:s}'.format( sampleNamesAr[i], sampleNamesAr[i+1] )
		# get the dictionaries
		if mDict1 == None:
			mFile1 = FileAllC_full( os.path.normpath( '{:s}/allc_{:s}{:s}.tsv'.format( allcPath, sampleNamesAr[i], ('_cov'+str(minCov) if minCov != None else '') ) ) )
			mDict1 = mFile1.getAllCDict( mtypes = methType )
			del(mFile1)
		
		mFile2 = FileAllC_full( os.path.normpath( '{:s}/allc_{:s}{:s}.tsv'.format( allcPath, sampleNamesAr[i+1], ('_cov'+str(minCov) if minCov != None else '') ) ) )
		mDict2 = mFile2.getAllCDict( mtypes = methType )
		if mDict1 == False or mDict2 == False:
			exit()
		#DMR analysis by chrm using processes
		#print(' analyzing by chrm with', numProc, 'processors' )
		results = [ pool.apply_async(processChrm, args=(mDict1[x], mDict2[x], dmrDict[x], comparisonName, x) ) for x in sorted(dmrDict.keys()) ]
		outArs = [p.get() for p in results]
		for a in outArs:
			outMat += a
		del(mDict1,mFile2)
		mDict1 = mDict2
	# end for
	del( mDict2 )
	# write output
	outFileStr = '{:s}_dmr-gen_{:s}.tsv'.format( outID, methType.lower() )
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outMat, info )
	print( 'Done' )

def readDMRFile( dmrFileStr ):
	'''
		return dictionary of subset (DMR) regions
		{chr1:[(start,end),(start2,end2)],chr2:[(start,end)],...}
	'''
	dmrFile = open( dmrFileStr, 'r' )
	dmrDict = {}
	dmrCount = 0
	
	for line in dmrFile:
			lineAr = line.rstrip().split()
			chrm = lineAr[0]
			start = int( lineAr[1] )
			end = int( lineAr[2] )
			dmrLabel = '{:s}:{:d}-{:d}'.format( chrm, start, end )
			if dmrDict.get(chrm) == None:
				dmrDict[chrm] = []
			dmrDict[chrm] += [(start, end, dmrCount, dmrLabel)]
			dmrCount += 1
	dmrFile.close()
	return dmrDict, dmrCount

def processChrm( mDict1, mDict2, dmrs, label, chrm ):
	print( ' analyzing', chrm )
	outAr = []
	for i in range(len(dmrs)):
		rStart, rEnd, rInt, rLabel = dmrs[i]
		#print( dmrs[i] )
		# analyze region
		regionStats = processRegion( rStart, rEnd, mDict1, mDict2 )
		rStatsAr = [ str(x) for x in regionStats ]
		tmpStr = '{:d}\t{:s}\t{:s}\t{:s}'.format( rInt, rLabel, label, '\t'.join( rStatsAr ) )
		outAr += [ tmpStr ]
	# end for
	#print( outAr )
	return outAr

def processRegion( start, end, mDict1, mDict2 ):
	cCount1, mCount1, tCount1 = regionCs( start, end, mDict1 )
	cCount2, mCount2, tCount2 = regionCs( start, end, mDict2 )
	cCount = min(cCount1, cCount2)
	tLength = end - start + 1
	# create output array
	wm1 = ( 0.0 if tCount1 == 0 else float(mCount1) / float(tCount1) )
	wm2 = ( 0.0 if tCount2 == 0 else float(mCount2) / float(tCount2) )
	return [ cCount, tLength, mCount1, mCount2, tCount1, tCount2, wm1, wm2, wm1-wm2]
		
def regionCs( start, end, inDict ):
	mCount = 0
	tCount = 0
	cCount = 0
	#print(len(inDict.keys()))
	for pos in range(start, end+1 ):
		tup = inDict.get( pos )
		if tup != None:
			cCount += 1
			mCount += tup[0]
			tCount += tup[1]
	# end for
	return cCount, mCount, tCount

'''def posCs( cSet, inDict ):
	mCount = 0
	tCount = 0
	for p in cSet:
		m,t = inDict[p]
		mCount += m
		tCount += t
	return mCount, tCount
'''
				
def writeOutput( outFileStr, outMat, info ):
	#print( outMat )
	header = info + '\nDMR\tregion\tlabel\tnum.cs\tlength\tmC.r1\tmC.r2\tt.r1\tt.r2\twm1\twm2\twm.diff\n'
	outFile = open( outFileStr, 'w' )
	outFile.write( header )
	for line in outMat:
		outFile.write( line + '\n' )
	outFile.close()
	
def parseInputs( argv ):
	numProc = 1
	outID = 'out'
	methType = 'C'
	minCov = COV
	startInd = 0
	
	for i in range(min(4,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			m = argv[i][3:].split(',')
			if len(m) > 1:
				print( 'WARNING: specify ONE methylation type at a time...choosing first listed' )
			methType = m[0]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum coverage must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	dmrFileStr = argv[startInd]
	allcPath = argv[startInd+1]
	if os.path.isdir( allcPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allcPath ) )
			exit()
	sampleNamesAr = []
	for i in range(startInd+2, len(argv)):
			sampleNamesAr += [ argv[i] ]
	if len( sampleNamesAr ) < 2:
		print( 'ERROR: must specify at least 2 generations' )
		exit()
	processInputs( dmrFileStr, allcPath, sampleNamesAr, outID, methType, numProc, minCov )
	

def printHelp():
	print( 'Usage: python dmr_gen_counts.py [-o=outID] [-m=methType] [-p=numProc] <dmrFile> <allcPath> <sample1> <sample2> [sampleN]*' )
	
if __name__ == "__main__":
	if len(sys.argv) < 5 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
