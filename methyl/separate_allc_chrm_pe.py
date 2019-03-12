import sys, math, glob, multiprocessing, subprocess, os, bisect, gzip

# Usage: python separate_allc_chrm_pe.py [-q] [-z] [-p=num_proc] <in_file> [in_file]*
# separate allc file into indiv files per chromosome
# output file names are derived from input file names

NUMPROC=1

def processInputs( inFileAr, numProc, isCompress, isPrint):
	
	if len(inFileAr) < numProc:
		numProc = len(inFileAr)
	if isPrint:
		print( 'Begin processing files with {:d} processors'.format( numProc ) )
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(allcFileStr, isCompress, isPrint) ) for allcFileStr in inFileAr ]
	
	suc = [ p.get() for p in results ]
	print( 'Done.' )

def processFile( allcFileStr, isCompress, isPrint ):
	
	if os.path.isfile( os.path.normpath( allcFileStr ) ) == False:
		return
	
	basename = os.path.basename( allcFileStr ).replace( '.gz', '' )
	info = '#from_script: separate_allc_chrm_pe.py; in_file: {:s}\n'.format( basename )
	
	outFileDict = {}
	curChrm = None
	curFile = None
	
	if isPrint:
		print( 'Processing', basename )
	
	if allcFileStr.endswith('.gz'):
		inFile = gzip.open(allcFileStr, 'rt' )
	else:
		inFile = open( allcFileStr, 'r' )
	
	for line in inFile:
		lineAr = line.rstrip().split('\t')
		if len(lineAr) < 6 or lineAr[5].isdigit() == False:
			continue
		chrm = lineAr[0]
		
		fName = outFileDict.get( chrm )
		# if this is a new chrm, close last file
		if chrm != curChrm and curFile != None:
			curFile.close()
			
		# need to create new file
		if fName == None:
			if isPrint:
				print( '  {:s} - {:s}'.format( basename, chrm ) )
			fName = getOutFileName( basename, chrm, isCompress )
			outFileDict[chrm] = fName
			if isCompress:
				curFile = gzip.open( fName, 'wt' )
			else:
				curFile = open( fName, 'w' )
		# need to reopen a file
		elif fName != None and chrm != curChrm:
			if isPrint:
				print( '  reopen {:s}'.format( fName ) )
			if isCompress:
				curFile = gzip.open(fName, 'at')
			else:
				curFile = open(fName, 'a')
		# write line
		curFile.write( line )
		curChrm = chrm
	# end for line
	
	inFile.close()
	# close output files
	#for key in outFileDict.keys():
	#	outFileDict[key].close()
	if curFile != None:
		curFile.close()
	
	if isPrint:
		print( 'Finished processing', basename )
	return True

def getOutFileName( basename, chrm, isCompress ):
	
	fInd = basename.rfind( '.' )
	if fInd == -1:
		tmp = basename + '_' + chrm + '.tsv'
	else:
		tmp = basename[:fInd] + '_' + chrm + basename[fInd:]
	return tmp + ('.gz' if isCompress else '')

def parseInputs( argv ):
	numProc = NUMPROC
	isCompress = False
	isPrint = True
	startInd = 0
	#print(min(4,len(argv)), len(argv))
	
	for i in range(min(4,len(argv))):
		if argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	inFileAr = []
	for j in range(startInd, len(argv)):
		inFileAr += [ argv[j] ]
	processInputs( inFileAr, numProc, isCompress, isPrint )

def printHelp():
	print ("\nUsage: python separate_allc_chrm_pe.py [-q] [-z] [-p=num_proc] <in_file> [in_file]*")
	print()
	print( 'Required:' )
	print( 'in_file\t\ttab-separated allc file(s)' )
	print()
	print( 'Optional:' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-z\t\tcompress allc output with gzip [default False]')
	print( '-p=num_proc\tnumber of processors to use [default 1]' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
