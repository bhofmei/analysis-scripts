import sys, math, multiprocessing, subprocess, os
from bioFiles import *

# Usage: python region_counts_pe.py [-h] [-q] [-k] [-o=outID] [-m=methTypes] [-p=numProc] [-v=minCov] <regionFile> <allcPath>  <sample1> <sample2> [sampleN]*
# compute weighted methylation for a set of regions
COV=None
NUMPROC=1
METHTYPE='C'
USEPICKLE=False

def processInputs( dmrFileStr, allcPath, sampleNamesAr, outID, methType, numProc, minCov, usePickle, isPrint ):
	if isPrint:
		print( 'Region file:', os.path.basename( dmrFileStr ) )
		print( 'AllC path:', allcPath )
		print( 'Samples:', ', '.join(sampleNamesAr) )
		print( 'Methylation type:',  methType )
		print( 'Minimum coverage:', minCov )
	info = '#from_script: region_counts_pe.py; dmr_file: {:s}; samples: {:s}; methyl_type: {:s}; min_cov: {:d}'.format( os.path.basename( dmrFileStr ), ','.join( sampleNamesAr ), methType, (minCov if minCov != None else 0) )

	# read DMR file
	if isPrint:
		print( 'Reading region file' )
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	chrmList = list(sorted(dmrDict.keys()))
	
	if isPrint:
		print( ' analyzing', dmrCount, 'regions' )

	outMat = []
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(handleChrm, args=(allcPath, sampleNamesAr, chrm, dmrDict[chrm], methType, minCov, isPrint, usePickle )) for chrm in chrmList ]
	outArs = [p.get() for p in results]
	for outAr in outArs:
		outMat += outAr
		
	# write output
	outFileStr = '{:s}_counts_{:s}.tsv'.format( outID, methType.lower() )
	if isPrint:
		print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outMat, info )
	if isPrint:
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
			dmrRegion = '{:s}:{:d}-{:d}'.format( chrm, start, end )
			dmrLabel = lineAr[3]
			if dmrDict.get(chrm) == None:
				dmrDict[chrm] = []
			dmrDict[chrm] += [(start, end, dmrCount, dmrRegion, dmrLabel)]
			dmrCount += 1
	dmrFile.close()
	return dmrDict, dmrCount

def handleChrm( allcPath, sampleNamesAr, chrm, dmrs, methType, minCov, isPrint, usePickle ):
	if isPrint:
		print( 'Processing', chrm )
	mDict1 = None
	outMat = []
	
	for i in range(len(sampleNamesAr)):
		smpName1 = sampleNamesAr[i]
		# get dictionaries
		mDict1 = getAllc(allcPath, smpName1, chrm, methType, minCov, isPrint, usePickle )
		
		# analyze dmrs
		chrmAr = processChrm( mDict1, dmrs, smpName1 )
		outMat += chrmAr
		
	# end for i
	return outMat

def getAllc( allcPath, sampleName, chrm, methType, minCov, isPrint, usePickle ):
	allcFileStr = os.path.normpath( '{:s}/allc_{:s}_{:s}{:s}.tsv'.format( allcPath, sampleName, chrm, ('_cov'+str(minCov) if minCov != None else '') ) )
	mFile = FileAllC_chrm( allcFileStr, isPrint=isPrint )
	if mFile.fileExists():
		mDict = mFile.getAllCDict( mtypes = [methType], isPickle=usePickle )
		return mDict
	else:
		return {}

def processChrm( mDict1, dmrs, label ):
	outAr = []
	for i in range(len(dmrs)):
		rStart, rEnd, rInt, rRegion, rLabel = dmrs[i]
		#print( dmrs[i] )
		# analyze region
		regionStats = processRegion( rStart, rEnd, mDict1 )
		rStatsAr = [ str(x) for x in regionStats ]
		tmpStr = '{:d}\t{:s}\t{:s}\t{:s}\t{:s}'.format( rInt, label, rRegion, rLabel, '\t'.join( rStatsAr ) )
		outAr += [ tmpStr ]
	# end for
	return outAr

def processRegion( start, end, mDict1 ):
	cCount, mCount, tCount = regionCs( start, end, mDict1 )
	tLength = end - start + 1
	# create output array
	wm = ( 0.0 if tCount == 0 else float(mCount) / float(tCount) )
	return [ cCount, tLength, mCount, tCount, wm]

def regionCs( start, end, inDict ):
	mCount = 0
	tCount = 0
	cCount = 0
	if inDict == {}:
		return cCount, mCount, tCount
	for pos in range(start, end+1 ):
		tup = inDict.get( pos )
		if tup != None:
			cCount += 1
			mCount += tup[0]
			tCount += tup[1]
	# end for
	return cCount, mCount, tCount

def writeOutput( outFileStr, outMat, info ):
	#print( outMat )
	header = info + '\nDMR\tsample\tregion\tlabel\tnum.cs\tlength\tmC.r\tt.r\twei_meth\n'
	outFile = open( outFileStr, 'w' )
	outFile.write( header )
	for line in outMat:
		outFile.write( line + '\n' )
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	outID = 'out'
	methType = METHTYPE
	minCov = COV
	isPrint = True
	usePickle = USEPICKLE
	startInd = 0

	for i in range(min(6,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-k':
			usePickle = True
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			m = argv[i][3:].split(',')
			if len(m) > 1:
				print( 'WARNING: specify ONE methylation type at a time...choosing first listed' )
			methType = m[0].upper()
			if methType not in [ 'CG', 'CHG', 'CHH', 'C' ]:
				print( 'ERROR: incorrect methylation type...choose one of "C", "CG", "CHG", "CHH"')
				exit()
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'WARNING: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
			startInd += 1
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
	processInputs( dmrFileStr, allcPath, sampleNamesAr, outID, methType, numProc, minCov, usePickle, isPrint )


def printHelp():
	print( 'Usage:\tpython region_counts_pe.py [-h] [-q] [-k] [-o=out_id]\n\t[-m=meth_type] [-p=num_proc] [-v=min_cov] <region_file> <allc_path> <sample1>\n\t<sample2> [sampleN]*' )
	print()
	print( 'Required:' )
	print( 'region_file\ttab-delimited file of DMR regions (1-indexed)' )
	print( 'allc_path\tpath to allC files' )
	print( 'sampleN\t\tsample name to analyze; order of samples is important' )
	print()
	print( 'Optional:' )
	print( '-h\t\tprint help and exit' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-k\t\tpickle allC files; creates/reads pickle version of allC files\n\t\tsaves subsequent computational time but is additional memory' )
	print( '-o=out_id\tidentifier for output file [default "out"]' )
	print( '-p=num_proc\tnumber of processors [default {:d}]'.format( NUMPROC) )
	print( '-m=meth_type\tsequence context; must be one of "CG", "CHG", "CHH", and "C"\n\t\t[default "C"]' )

if __name__ == "__main__":
	if len(sys.argv) < 5 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
