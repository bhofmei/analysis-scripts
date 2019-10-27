import sys, math, glob, multiprocessing, subprocess, os, bisect, gzip
import pysam #, pybedtools

NBINS=10
BINSIZE=100
NUMPROC=1

# Usage: python heatmap_abs_coverage_pe.py [-h] [-q] [-z] [-b=num_bins | -s=bin_size] [-o=out_id] [-l=labels] [-p=num_proc] [-m=min_length] <region_file> <reads_file> [reads_file]*

def processInputs( regionFileStr, readFileAr, numBins, binSize, outId, labels, minLength, maxLength, numProc, isPrint, isCompress, isCentered, isScaled ):
	
	# get sample names
	if labels == None:
		labels = [ getSampleName(file) for file in readFileAr ]
		
	if isPrint:
		print('Region file:', os.path.basename(regionFileStr))
		print( 'Read files: ', ', '.join([os.path.basename(file) for file in readFileAr]) )
		print( 'Min region length:', minLength )
		print( 'Max region length:', maxLength )
	regionAr, kLength = readRegions( regionFileStr, minLength, maxLength )
	#print(kLength)
	# determine bin number and bin size
	if binSize == -1:
		binSize = int( math.ceil(float(kLength) / numBins) )
	else:
		numBins = int( math.ceil(float(kLength) / binSize) )
	
	if isPrint:
		print( 'Number of regions:', len(regionAr) )
		print( numBins, 'bins of size', binSize )
	
	info = 'from_script: heatmap_abs_coverage_pe.py; region_file: '
	info += os.path.basename(regionFileStr) + '; bin_size: ' + str(binSize)
	info += '; num_bins: ' + str(numBins) + '; is_scaled: ' + str(isScaled)
	info += '; min_size: ' + str(minLength) + '; max_length: ' + str(maxLength)
	
	if isPrint:
		print( 'Begin processing {:d} samples with {:d} processors'.format(len(labels), numProc) )
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(processFile, args=(readFileAr[i], regionAr, labels[i], binSize, numBins, outId, info, isPrint, isCompress, isCentered, isScaled) ) for i in range(len(readFileAr)) ]
	out = [p.get() for p in results]
	print('Done')
		
def getSampleName( fileStr ):
	bName = os.path.basename(fileStr)
	rIndex = bName.rfind('.')
	if rIndex != -1:
		return bName[:rIndex]
	else:
		return bName

def readRegions( regionFileStr, minLength, xLength ):
	outAr = []
	
	bedFile = open( regionFileStr, 'r' )
	maxLength = -1
	count = 1
	
	for line in bedFile:
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) start (2) end (3) name?
		if len(lineAr) < 3:
			continue
		chrm = lineAr[0]
		start = int(lineAr[1])
		end = int(lineAr[2])
		regionLength = end - start
		if regionLength >= minLength and regionLength <= xLength:
			name = ( 'region-'+count if len(lineAr) < 4 else lineAr[3] )
			outAr += [(chrm, start, end, name)]
			count += 1
			maxLength = max(regionLength, maxLength)
	# end for line
	bedFile.close()
	return outAr, maxLength
	
def processFile( fileStr, regionAr, label, binSize, numBins, outId, info, isPrint, isCompress, isCentered, isScaled ):
	
	# file is a bam file
	if fileStr.endswith('.bam'):
		outMatrix = processBAM( fileStr, regionAr, binSize, numBins, isCentered, isScaled )
	else:
		outMatrix = processBED()
	
	outFileStr = outId + '_' if outId != None else ''
	outFileStr += label + '_abs_heatmap.tsv'
	outFileStr += '.gz' if isCompress else ''
	
	if isPrint:
		print('Output of {:s} written to {:s}'.format( os.path.basename(fileStr), outFileStr) )
	
	writeOutput(outFileStr, outMatrix, regionAr, info, isCompress)

def checkBAM( bamFileStr ):
	bamIndex = bamFileStr + '.bai'
	if os.path.isfile(bamFileStr) == False:
		print('ERROR: BAM file does not exist')
		return False
	elif os.path.isfile(bamIndex) == False:
		print('WARNING: BAM index file does not exist...creating')
		pysam.index( bamFileStr )
	return True	

def processBAM( fileStr, regionAr, binSize, numBins, isCentered, isScaled ):
	# check file exists and bam index
	isValid = checkBAM( fileStr )
	if not isValid:
		return
	
	# otherwise open bam file
	bamFile = pysam.AlignmentFile( fileStr, 'rb' )
	scaleVal = float(bamFile.mapped / 1000000) if isScaled else 1.0
	
	outMatrix = []
	
	for region in regionAr:
		covAr = processBAMRegion( bamFile, region, binSize, numBins, isCentered, scaleVal )
		outMatrix.append(covAr)
	# end for region
	bamFile.close()
	return outMatrix

def processBAMRegion( bamFile, region, binSize, numBins, isCentered, scaleVal ):
	
	outAr = [-1] * numBins
	chrm, start, end, label = region
	
	if isCentered:
		midPoint = int( (end - start)/2.0 ) + start
		lBins = int( math.ceil((midPoint - start) / binSize ))
		actBins = lBins * 2
		startBin = int (numBins / 2) - lBins
	# end if isCentered
	else:
		actBins = int( math.ceil((end - start)/binSize) )
		startBin = 0
	
	for i in range(actBins):
		pStart = i*binSize + start
		pEnd = pStart + binSize
		pBin = startBin + i
		regionVal = computeBAMRegion( bamFile, chrm, pStart, pEnd )
		outAr[pBin] = regionVal / scaleVal
	# end for i
	return outAr
			
def computeBAMRegion( bamFile, chrm, start, end ):
	reads = bamFile.fetch(chrm, start, end )
	readsAr = [1 for read in reads]
	return sum(readsAr)

def writeOutput( outFileStr, outMatrix, regionAr, info, isCompress ):
	
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	headerAr = ['chrm', 'start', 'end', 'region', 'length', 'bin', 'value']
	outFile.write( info + '\n' + '\t'.join(headerAr) + '\n')
	
	for i in range(len(outMatrix)):
		chrm, start, end, name = regionAr[i]
		regionStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:d}'.format( chrm, start, end, name, end-start )
		tmpAr = outMatrix[i]
		for j in range(len(tmpAr)):
			valStr = '\t{:d}\t'.format(j)
			valStr += 'NA' if tmpAr[j] == -1 else '{:.4f}'.format(tmpAr[j])
			outFile.write( regionStr + valStr + '\n')
	# end for x
	outFile.close()
		
def parseInputs( argv ):
	isPrint = True
	isCompress = False
	isCentered = False
	isScaled = True
	numBins = -1
	binSize = -1
	minLength = 0
	maxLength = math.inf
	outId = None
	labels = None
	numProc = NUMPROC
	startInd = 0
	
	for i in range(len(argv)):
		if argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i] == '-c':
			isCentered = True
			startInd += 1
		elif argv[i] == '-n':
			isScaled = False
			startInd += 1
		elif argv[i].startswith( '-b=' ):
			# check for previous '-s'
			if binSize != -1:
				print( 'ERROR: cannot specify -b and -s together' )
				exit()
			try:
				numBins = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of bins must be an integer' )
				exit()
		elif argv[i].startswith( '-s=' ):
			# check previous '-b'
			if numBins != -1:
				print( 'ERROR: cannot specify -b and -s together' )
				exit()
			try:
				binStr = argv[i][3:]
				if binStr.endswith( 'k' ) or binStr.endswith( 'K' ):
					binSize = int( binStr[:-1] ) * 1000
				elif binStr.endswith( 'm' ) or binStr.endswith( 'M' ):
					binSize = int( binStr[:-1] ) * 1000000
				else:
					binSize = int( binStr )
				startInd += 1
			except ValueError:
				print( 'ERROR: bin size must be an integer' )
				exit()
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i].startswith( '-m=' ):
			try:
				minLength = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: min length must be an integer' )
				exit()
		elif argv[i].startswith( '-x=' ):
			try:
				maxLength = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: max length must be an integer' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outId = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-l='):
			labels = argv[3:].split(',')
			startInd += 1
		elif argv[i].startswith('-'):
			print( 'ERROR: unknown parameter {:s}. Use -h to see available parameters'.format( argv[i] ) )
			exit()
	# end for i
	
	if numBins == -1 and binSize == -1:
		binSize = BINSIZE
	
	regionFileStr = argv[startInd]
	readFileAr = []
	
	for i in range(startInd + 1, len(argv)):
		readFileAr += [ argv[i] ]
	
	if labels != None and len(readFileAr) != len(labels):
		print('WARNING: number of labels does not match number of read files. Ignoring labels parameter')
		labels = None
	
	processInputs( regionFileStr, readFileAr, numBins, binSize, outId, labels, minLength, maxLength, numProc, isPrint, isCompress, isCentered, isScaled )

def printHelp():
	print()
	print('Usage:\tpython heatmap_abs_coverage_pe.py [-h] [-q] [-z] [-c] [-n]')
	print('\t[-b=num_bins | -s=bin_size] [-o=out_id] [-l=labels] [-m=min_length]')
	print('\t[-x=max_length] [-p=num_proc] <region_file> <reads_file> [reads_file]*')
	print()
	print('Required:')
	print('region_file\tBED file with regions for analysis')
	print('reads_file\tBED or BAM file with mapped reads')
	print('\t\tmultiple reads files can be specified as individual samples')
	print()
	print('Optional:')
	print('-h\t\tprint this help message and exit')
	print('-q\t\tquiet; do not print progress')
	print('-z\t\tcompress output with gzip')
	print('-c\t\tcenter the heatmap output [default left-justified]')
	print('-n\t\tdo NOT scale coverage counts by library size [default is scaled]')
	print('-b=num_bins\tnumber of bins for map; size of bin determined by longest region')
	print('-s=bin_size\tsize of bin; number of bins deteremined by longest region\n\t\t[default {:d} bp]'.format(BINSIZE))
	print('-o=out_id\toutput file prefix [default read file name]')
	print('-l=labels\tcomma-separated labels for samples to use in output file name\n\t\t[default read file name]')
	print('-m=min_length\tminimum length of regions to include in output\n\t\t[default include all]')
	print('-x=max_length\tmaximum length of regions to include in output\n\t\t[default include all]')
	print('-p=num_proc\tnumber of processors [default 1]')
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
