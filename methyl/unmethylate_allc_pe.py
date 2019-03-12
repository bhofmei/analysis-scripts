import sys, multiprocessing, subprocess, os, gzip

# Usage: python unmethylate_allc_pe.py [-f] [-q] [-h] [-z] [-p=num_proc] [-v=coverage] <allc_file> [allc_file]*
# creates a pseudo-allc file where "reads" were unmethylated
# uses the same per-position coverage as input file unless otherwise specified

NUMPROC=1

def processInputs( allCFileAr, numProc, newCov, sampleFile, isCompress, isPrint ):

	if sampleFile:
		if isPrint:
			print( 'Reading file of samples names...' )
		allCFileAr = readSampleFile( allCFileAr[0] )

	tmpAr = [ os.path.basename(x) for x in allCFileAr ]
	if isPrint:
		print( 'Files: {:s}\nCoverage: {:s}'.format( ' '.join( tmpAr ), ( 'as-is' if newCov == None else str(newCov) ) ) )

	# check for all allC files before doing any work
	for allcFile in allCFileAr:
		if checkFiles( allcFile ) == False:
			exit()

	#  loop through files
	if isPrint:
		print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(allcFile, newCov, isCompress, isPrint) ) for allcFile in allCFileAr ]
	suc = [ p.get() for p in results ]
	if isPrint:
		print( 'Done' )

def readSampleFile( fileStr ):

	if os.path.isfile(fileStr) == False:
		print( 'ERROR:', os.path.basename(fileStr), 'does not exist' )
		exit()
	sampleAr = []
	inFile = open( fileStr, 'r' )
	for line in inFile:
		name = line.rstrip()
		sampleAr += [ name ]
	inFile.close()
	return sampleAr

def checkFiles( allcFile ):
	if os.path.isfile( allcFile ) == False and os.path.isfile( allcFile +'.gz' ) == False:
		print( 'ERROR: allC file {:s} not found'.format( os.path.basename( allcFile) ) )
		return False
	return True

def processFile( allcFileStr, newCov, isCompress, isPrint ):
	outFileStr = getOutFileName( allcFileStr, isCompress )
	if isPrint:
		print( '  {:s} -> {:s}'.format( os.path.basename(allcFileStr), os.path.basename(outFileStr) ) )

	allCFile = openFile( allcFileStr )
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	# read file
	for line in allCFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		if newCov != None:
			lineAr[5] = str(newCov)
		# replace mc_count and methylated
		lineAr[4] = '0'
		lineAr[6] = '0'
		outFile.write( '\t'.join( lineAr ) + '\n' )
	# end for line
	allCFile.close()
	outFile.close()
	return True

def getOutFileName( allcFileStr, isCompress ):
	bName = os.path.basename( allcFileStr )
	bName = bName.replace('.gz','')
	rInd = bName.rfind( '_' )
	if rInd == -1:
		cName = bName.replace('.tsv','-unmethylated.tsv')
	else:
		cName = bName[:rInd] + '-unmethylated' + bName[rInd:]
	cName += '.gz' if isCompress else ''
	return os.path.join( os.path.dirname(allcFileStr), cName )


def openFile( allcFileStr ):
	if os.path.isfile(allcFileStr) and allcFileStr.endswith('.gz'):
		return gzip.open( allcFileStr, 'rt' )
	elif os.path.isfile( allcFileStr ) == False and os.path.isfile( allcFileStr + '.gz' ):
		return gzip.open( allcFileStr+'.gz', 'rt' )
	else:
		return open( allcFileStr, 'r' )

def parseInputs( argv ):
	numProc = NUMPROC
	sampleFile = False
	newCov = None
	isPrint = True
	isCompress = False
	startInd = 0

	for i in range(min(6,len(argv)-1)):
		if argv[i] == '-f':
			sampleFile = True
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
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
				newCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: new position coverage must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for

	allcFilesAr = []
	for j in range( startInd, len(argv) ):
		allcFilesAr += [ argv[j] ]
		if j == 1 and sampleFile:
			print( 'ERROR: only specify one file with sample names' )
			exit()

	processInputs( allcFilesAr, numProc, newCov, sampleFile, isCompress, isPrint )

def printHelp():

	print( 'Usage:\tpython unmethylate_allc_pe.py [-q] [-h] [-f] [-p=num_proc] [-v=coverage]\n\t<allc_file> [allc_file]*' )
	print( )
	print( 'Required:' )
	print( 'allc_file\tallC file to unmethylated\n\t\twhen "-f" set, file with list of allC files' )
	print()
	print( 'Optional:' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-h\t\tprint help message and exit' )
	print( '-z\t\tcompress output with gzip')
	print( '-f \t\tallC files names listed in the file' )
	print( '-p=num_proc\tnumber of processors to use [default {:d}]'.format( NUMPROC) )
	print( '-v=coverage\tcoverage for each position [default as-is in input]' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
