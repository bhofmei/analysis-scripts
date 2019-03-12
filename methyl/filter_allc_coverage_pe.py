import sys, math, glob, multiprocessing, subprocess, os, gzip
from bioFiles import *

# Usage: python filter_allc_coverage_pe.py [-z] [-v=min_cov] [-num_proc] [-c=chrm_list | -cf=fasta_index] <allc_path> <sample1> [sampleN]*
# creates new allc files for all input samples that only includes information
# about positions which have at least minCov reads for each sample
# expects all chromosomes in one allC file

PICKLE=False
MINCOV=3
NUMPROC=1
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def processInputs( allcPath, sampleNamesAr, chrmList, fastaIndex, minCov, numProc, isCompress, isPrint ):
	if isPrint:
		print( 'AllC path:', allcPath )
		print( 'Samples:', ', '.join(sampleNamesAr) )
		print( 'Minimum coverage:', minCov )

	info = '#from_script: filter_allc_coverage_pe.py; min_cov: {:d}; samples_included: {:s}\n'.format( minCov, ','.join(sampleNamesAr) )

	if fastaIndex != None:
		chrmList = readFastaIndex( fastaIndex )
	
	# loop through all samples and get sets of positions with minCov
	if isPrint:
		print( 'Analyzing {:d} chromosomes of {:d} samples with {:d} processes'.format( len(chrmList), len(sampleNamesAr), numProc) )
		
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(processChrm, args=(allcPath, sampleNamesAr, chrm, minCov, info, isCompress, isPrint) ) for chrm in chrmList ]
	pResults = [p.get() for p in results]
	print( 'Done' )

def readFastaIndex( fastaIndexStr ):
	if os.path.isfile(fastaIndexStr) == False:
		print( 'ERROR:', os.path.basename(fastaIndexStr), 'does not exist' )
		exit()
	chrmList = []
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		chrmList += [ lineAr[0] ]
	fastaIndex.close()
	return chrmList

def processChrm( allcPath, sampleNamesAr, chrm, minCov, info, isCompress, isPrint):
	if isPrint:
		print( 'Processing chrm {:s} for {:d} samples'.format(chrm, len(sampleNamesAr)))
	# get the positions
	posSet = getCovPositions( allcPath, sampleNamesAr, chrm, minCov, isPrint)
	
	if len(posSet) == 0:
		print( 'WARNING: no positions available for chrm', chrm )
		return False
	else:
		# write to files
		for sampleName in sampleNamesAr:
			writeAllcFile( allcPath, sampleName, chrm, posSet, minCov, info, isCompress, isPrint )
		# end for sampleName
		return True

def getCovPositions( allcPath, sampleNamesAr, chrm, minCov, isPrint ):
	# loop through samples
	outSet = []
	for i in range(len(sampleNamesAr)):
		sampleName = sampleNamesAr[i]
		allcFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sampleName, chrm))
		posSet = getCovPositionSample( allcFileStr, minCov, isPrint )
		# first sample starts the set
		if i == 0:
			outSet = posSet
		else:
			# do intersection
			outSet = outSet & posSet
		# end if i == 0
		if isPrint:
			print('-sample {:s}, positions: {:d}, intersection: {:d}'.format(sampleName, len(posSet), len(outSet)))
	# end for i
	return outSet

def getCovPositionSample( allcFileStr, minCov, isPrint ):
	posAr = []
	mFile = FileAllC_chrm( allcFileStr, isPrint=isPrint )

	if mFile.fileExists():
		mDict = mFile.getAllCDict( mtypes = ['C'], isPickle=PICKLE)
		for pos in mDict.keys():
			# check coverage
			try:
				if mDict[pos][1] >= minCov:
					posAr += [ pos ]
			except KeyError:
				print( 'KEY ERROR: ', pos)
		# end for pos
	# end if
	return set( posAr )

def combineCovSamples( posDictAr ):
	# posDictAr is array of dictionaries with chrms
	chrmList = list( posDictAr[0].keys() )
	outDict = {}
	# loop through chrm
	for chrm in chrmList:
		# start with set 0 and intersection remaining
		cSet = posDictAr[0][chrm]
		for i in range(1, len(posDictAr) ):
			newSet = posDictAr[i][chrm]
			cSet = cSet & newSet
		# end for i
		outDict[chrm] = cSet
	return outDict

def writeAllcFile( allcPath, sampleName, chrm, posSet, minCov, info, isCompress, isPrint ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sampleName, chrm ) )
	inFileObj = FileBio( inFileStr )
	inFile = inFileObj.fbOpen()
	outFileStr = inFileObj.fbBasename() + '_cov{:d}.tsv'.format( minCov )
	if isCompress:
		outFileStr += '.gz'
	if isPrint:
		print( 'Writing output to', os.path.basename( outFileStr ) )
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write( info )

	for line in inFile:
		if line.startswith( '#' ):
			outFile.write( line )
			continue
		lineAr = line.rstrip().split( '\t' )
		if len( lineAr ) < 7 or lineAr[6].isdigit() == False:
			continue
		pos = int( lineAr[1] )
		if pos in posSet:
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()

def parseInputs( argv ):
	minCov = MINCOV
	numProc = NUMPROC
	isCompress = False
	chrmList = None
	fastaIndex = None
	isPrint = True
	startInd = 0

	for i in range(min(6,len(argv))):
		if argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )

			except ValueError:
				print( 'WARNING: minimum coverage must be integer...using default', MINCOV )
				minCov = MINCOV
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'WARNING: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i].startswith( '-c='):
			if fastaIndex != None:
				print( 'ERROR: cannot specify chromosome list and fasta index file' )
				exit()
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-cf=' ):
			if chrmList != None:
				print( 'ERROR: cannot specify chromosome list and fasta index file' )
				exit()
			fastaIndex = argv[i][4:]
			startInd += 1	
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	allcPath = argv[startInd]
	if os.path.isdir( allcPath ) == False:
		print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allcPath ) )
		exit()
	sampleNamesAr = []
	for i in range(startInd+1, len(argv)):
		sampleNamesAr += [ argv[i] ]
	processInputs( allcPath, sampleNamesAr, chrmList, fastaIndex, minCov, numProc, isCompress, isPrint )

def printHelp():
	print( 'Usage:\tpython filter_allc_coverage_pe.py [-q] [-h] [-z] [-v=min_cov] [-c=chrm_list | -cf=fasta_index]\n\t<allc_path> <sample_name> [sample_name]*' )
	print()
	print( 'Required:' )
	print( 'allc_path\tpath to allC files' )
	print( 'sampleN\t\tname of sample; used to find allC files' )
	print()
	print( 'Optional:' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-h\t\tprint help and exit' )
	print( '-z\t\tcompress allc output with gzip [default False]')
	print( '-v=min_cov\tmin coverage for positions to include [default {:d}]'.format( MINCOV) )
	print( '-p=num_proc\tnumber of processors to use [default {:d}]'.format( NUMPROC) )
	print( '-c=chrm_list\tcomma-separated list of chrms to use [default for arabidopsis]' )
	print( '-cf=fasta_index\tfasta index file with chrms to use' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
