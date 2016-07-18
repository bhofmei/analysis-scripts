import sys, multiprocessing, subprocess, os, gzip

# Usage: python unmethylate_allc_pe.py [-f] [-p=num_proc] [-v=NA] <allc_file> [allc_file]*
# creates a pseudo-allc file where "reads" were unmethylated
# uses the same per-position coverage as input file unless otherwise specified

NUMPROC=1

def processInputs( allCFileAr, numProc, newCov, sampleFile ):
	
	if sampleFile:
		print( 'Reading file of samples names...' )
		allCFileAr = readSampleFile( allCFileAr[0] )
		
	tmpAr = [ os.path.basename(x) for x in allCFileAr ]
	print( 'Files: {:s}\nCoverage: {:s}'.format( ' '.join( tmpAr ), ( 'as-is' if newCov == None else str(newCov) ) ) )
	
	# check for all allC files before doing any work
	for allcFile in allCFileAr:
		if checkFiles( allcFile ) == False:
			exit()
	
	#  loop through files
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(allcFile, newCov) ) for allcFile in allCFileAr ]
	suc = [ p.get() for p in results ]
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

def processFile( allcFileStr, newCov ):
	outFileStr = getOutFileName( allcFileStr )
	print( '  {:s} -> {:s}'.format( os.path.basename(allcFileStr), os.path.basename(outFileStr) ) )
	
	allCFile = openFile( allcFileStr )
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
	
def getOutFileName( allcFileStr ):
	bName = os.path.basename( allcFileStr )
	#bName = bName.replace('allc_','')
	rInd = bName.rfind( '_' )
	if rInd == -1:
		cName = bName.replace('.tsv','-unmethylated.tsv')
	else:
		cName = bName[:rInd] + '-unmethylated' + bName[rInd:]
	return os.path.join( os.path.dirname(allcFileStr), cName )
	 

def openFile( allcFileStr ):
	if os.path.isfile( allcFileStr ) == False and os.path.isfile( allcFileStr + '.gz' ):
		return gzip.open( allcFileStr+'.gz', 'rt' )
	else:
		return open( allcFileStr, 'r' )

def parseInputs( argv ):
	numProc = NUMPROC
	sampleFile = False
	newCov = None
	startInd = 0
	
	for i in range(min(5,len(argv)-2)):
		if argv[i] == '-f':
			sampleFile = True
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
	
	processInputs( allcFilesAr, numProc, newCov, sampleFile )

def printHelp():
	
	print( 'Usage: python unmethylate_allc_pe.py [-f] [-p=num_proc] [-v=NA] <allc_file> [allc_file]*' )
	print( 'Creates pseudo-allC file with all positions unmethylated' )
	print( 'Any or all input allc files may be compressed with gzip (file name ends in .gz)' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
