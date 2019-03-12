import sys, math, glob, multiprocessing, subprocess, os, random

# Usage: python3.4 subset_coverage_from_bam.py [-keep] [-p=num_proc] [-o=out_prefix] <-g=genome_size | -f=fasta_index> <coverage_values> <bam_file | sam_file >

def processInputs( bamFileStr, coverageValues, genomeSize, fastaIndexStr, outPre, numProc, rmTmp ):
	pool = multiprocessing.Pool( processes=numProc )
	if genomeSize == -1:
		print( 'Reading fasta index...' )
		genomeSize = readFastaIndex( fastaIndexStr )
	if outPre == None:
		baseName = os.path.basename(bamFileStr)
		outPre = baseName[:(baseName.rfind('.'))]
	# possibly convert to sam
	print( 'Checking if {:s} is a bam file...'.format( bamFileStr) )
	tmpFileStr = checkBam( bamFileStr )
	# sam file
	if tmpFileStr == False:
		isBam = False
	else:
		bamFileStr = tmpFileStr
		isBam = True
	
	# check size
	print( 'Analyzing {:s}...'.format( bamFileStr ) )
	numReads, readLength = checkFileSize( bamFileStr )
	if float(numReads*readLength) < genomeSize*max(coverageValues):
		print( 'ERROR: not enough reads in bam/sam file\nTry a new file or a smaller coverage value' )
		exit()
	print( 'Generating {:d} random sets of number using {:d} processors'.format( len(coverageValues), numProc ) )
	results = [ pool.apply_async( generateRandomSet, args=(x, numReads, genomeSize, readLength)) for x in coverageValues ]
	lineSetsAr = [ p.get() for p in results ]
	lineNumsAr, outFileAr, outFileStrAr = prepReading( lineSetsAr, outPre, coverageValues )
	print( 'Reading {:s}...'.format( bamFileStr ) )
	readSam( bamFileStr, lineNumsAr, outFileAr )
	
	print( 'Converting subset sam files to bam with {:d} processors...'.format( numProc) )
	
	results = [ pool.apply_async( convertToBam, args=(f, rmTmp)) for f in outFileStrAr ]
	f = [ p.get() for p in results ]
	
	# remove temporary sam if necessary
	if isBam and rmTmp == False:
		os.remove( tmpFileStr )
	print( 'Done.' )

def readFastaIndex( fastaIndexStr ):
	
	genomeSize = 0
	excludeList = ['ChrC', 'ChrM', 'chloroplast', 'mitochondria' ]
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		# possible headers
		if line.startswith( '#' ) or line.startswith( '@' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) size
		if lineAr[0] not in excludeList:
			genomeSize += int( lineAr[1] )
	fastaIndex.close()
	return genomeSize

def checkFileSize( samFileStr ):
	command = "grep -vc '^@' {:s}".format( samFileStr )
	tmpReads = subprocess.check_output( command, shell=True, universal_newlines = True )
	command = "grep -v '^@' {:s} | head -n 10".format( samFileStr )
	lines = subprocess.check_output( command, shell=True, universal_newlines = True )
	sum = 0
	n = 0
	for line in lines.split('\n'):
		lineAr = line.rstrip().split('\t')
		if len( lineAr ) < 10:
			continue
		sum += len( lineAr[9] )
		n += 1
	rLen = float(sum) / n
	return int( tmpReads.split()[0] ), rLen

def generateRandomSet( coverageValue, numReads, genomeSize, readLength ):
	
	readsNeeded = coverageValue * genomeSize / readLength
	currSet = set()
	while len(currSet) < readsNeeded:
		# generate random read number
		t = random.randrange(0,numReads)
		currSet.add(t)
	return currSet

def prepReading( lineSetsAr, outPre, coverageValues ):
	n = len( coverageValues )
	lineNumsAr = []
	for s in lineSetsAr:
		lineNumsAr.append( sorted(list(s)) )
	outFileAr = []
	outFileStrAr = []
	for i in range(n):
		outFileStr = '{:s}_cov{:.1f}.sam'.format( outPre, coverageValues[i] )
		outFile = open( outFileStr, 'w' )
		outFileAr += [ outFile ]
		outFileStrAr += [ outFileStr ]
	return lineNumsAr, outFileAr, outFileStrAr

def checkBam( bamFileStr ):
	if bamFileStr.endswith( '.bam' ):
		print( 'Converting to sam...' )
		tmpStr = bamFileStr.replace( '.bam','.tmp.sam' )
		command = 'samtools view -h {:s} > {:s}'.format( bamFileStr, tmpStr )
		subprocess.call( command, shell = True )
		return tmpStr
	return False

def readSam( samFileStr, lineNumsAr, outFileAr ):
	
	samFile = open( samFileStr, 'r' )
	n = len( lineNumsAr )
	curInd = [0] * n
	curLine = [ lineNumsAr[i][0] for i in range(n)]
	count = 0
	
	# read file
	for line in samFile:
		# headers get written to all files
		if line.startswith( '@' ):
			for i in range(n):
				outFileAr[i].write( line )
			continue
		# current line not needed
		if count not in curLine:
			count += 1
		# at least one file needs this line
		else:
			# loop through curLine to determine file
			for i in range(n):
				if curLine[i] == count:
					# update
					outFileAr[i].write( line )
					curInd[i] += 1
					if curInd[i] < len( lineNumsAr[i] ):
						curLine[i] = lineNumsAr[i][curInd[i]]
					else:
						curLine[i] = -1
			# end for i
			count += 1
		if sum( curInd ) == -1*n:
			break
	# end for line
	samFile.close()
	# close out files
	y = [ x.close() for x in outFileAr ]
	return True

def convertToBam( samFileStr, removeTmp ):
	
	command = 'samtools view -c {:s}'.format( samFileStr )
	tmpReads = subprocess.check_output( command, shell=True, universal_newlines = True )
	print( '{:s} - {:s} reads'.format( samFileStr, tmpReads.split()[0] ) )
	bamFileStr = samFileStr.replace( '.sam','.bam' )
	command = 'samtools view -b -o {:s} {:s}'.format( bamFileStr, samFileStr )
	subprocess.call( command, shell=True )
	command = 'samtools index {:s}'.format( bamFileStr )
	subprocess.call( command, shell=True )
	
	if removeTmp:
		os.remove( samFileStr )
	
def parseInputs( argv ):
	
	outPre = None
	genomeSize = -1
	fastaIndexStr = None
	startInd = 0
	numProc = 1
	rmTmp = True
	
	for i in range(min(5,len(argv)-2)):
		if argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERRROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-f=' ):
			if genomeSize != -1:
				print('ERROR: cannot specify fasta index file and genome size')
				exit()
			fastaIndexStr = argv[i][3:]
			print( fastaIndexStr )
			startInd += 1
		elif argv[i] == '-keep':
			rmTmp = False
			startInd += 1
		elif argv[i].startswith( '-g=' ):
			if fastaIndexStr != None:
				print('ERROR: cannot specify fasta index file and genome size')
				exit()
			genomeSize = argv[i][3:]
			try:
				if genomeSize.endswith( 'k' ) or genomeSize.endswith( 'K' ):
					genomeSize = int( genomeSize[:-1] ) * 1000
				elif genomeSize.endswith( 'm' ) or genomeSize.endswith( 'M' ):
					genomeSize = int( genomeSize[:-1] ) * 1000000
				elif genomeSize.endswith( 'g' ) or genomeSize.endswith( 'G' ):
					genomeSize = int( genomeSize[:-1] ) * 1000000000
				else:
					genomeSize = int( genomeSize )
			except ValueError:
				print( 'Genome size must be integer' )
				exit()
			startInd += 1
	# genome size or fasta file must be specified
	if genomeSize == -1 and fastaIndexStr == None:
		print( 'ERROR: genome size or fasta index file must be specified' )
		exit()
	# end for
	coverage = argv[startInd].split( ',' )
	try:
		coverageValues = [ float(x) for x in coverage ]
	except ValueError:
		print( 'ERROR: at least one coverage value not numeric' )
		exit()
	bamFileStr = argv[startInd + 1]
	processInputs( bamFileStr, coverageValues, genomeSize, fastaIndexStr, outPre, numProc, rmTmp )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 subset_coverage_from_bam.py [-keep] [-p=num_proc] [-o=out_prefix] <-g=genome_size | -f=fasta_index> <coverage_values> <bam_file | sam_file >")
	else:
		parseInputs( sys.argv[1:] )
