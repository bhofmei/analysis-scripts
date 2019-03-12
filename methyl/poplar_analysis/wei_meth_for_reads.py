import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import pysam

WINDSIZE=100
MINCOV=0.5

# Usage: python3 wei_meth_for_reads.py [-h] [-q] [-o=outID] [-w=window_size] [-v=min_overlap] <fasta_file> <bam_file>

def processInputs( fastaFileStr, bamFileStr, outID, windowSize, minCov, isPrint ):
	if isPrint:
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'Window size:', windowSize )
		print( 'Min percent overlap:', minCov )
	
	# check BAM file/index
	checkBAM( bamFileStr )
	# check FASTA file/index
	checkFASTA( fastaFileStr )
	
	# get genome regions
	genome = getGenome( fastaFileStr )
	
	# create pysam file objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	fastaFile = pysam.Fastafile( fastaFileStr )
	
	# process
	if isPrint:
		print('Analyzing {:d} chrms'.format( len(genome)))
		
	outMatrix = []
	for i in range(len(genome)):
		chrmRes = processChrm( fastaFile, bamFile, genome[i], windowSize, minCov, isPrint )
		outMatrix.append(chrmRes)
	
	bamFile.close()
	fastaFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_mReads.tsv'
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: wei_meth_for_reads.py; window_size:{:d}; min_cov:{:g}\n'.format(windowSize, minCov)
	
	writeOutput( outFileStr, outMatrix, genome, info )
	print('Done')

def checkBAM( bamFileStr ):
	bamIndex = bamFileStr + '.bai'
	if os.path.isfile(bamFileStr) == False:
		print('ERROR: BAM file does not exist')
		exit()
	elif os.path.isfile(bamIndex) == False:
		print('WARNING: BAM index file does not exist...creating')
		pysam.index( bamFileStr )
	return True

def checkFASTA( fastaFileStr ):
	fastaIndex = fastaFileStr + '.fai'
	if os.path.isfile(fastaFileStr) == False:
		print('ERROR: FASTA file does not exist')
		exit()
	elif os.path.isfile(fastaIndex) == False:
		print('WARNING: FASTA index file does not exist...creating')
		pysam.faidx( fastaFileStr )
	return True

def getGenome( fastaFileStr ):
	fastaIndexStr = fastaFileStr + '.fai'
	
	outAr = []
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) length
		chrm = lineAr[0]
		chrmLen = int(lineAr[1])
		outAr += [(chrm, chrmLen)]
	# end for
	fastaIndex.close()
	
	return outAr
	
def processChrm( fastaFile, bamFile, chrmInfo, windowSize, minCov, isPrint ):
	
	if isPrint:
		print('Analyzing chrm', chrmInfo[0] )
		
	chrmName = chrmInfo[0]
	chrmLen = chrmInfo[1]
	nBins = int( math.floor(float(chrmLen)/ windowSize) )
	
	outAr = []
	# loop through regions
	for i in range(nBins-1):
		start = i * windowSize
		end = (i+1) * windowSize
		regionVals = processRegion( fastaFile, bamFile, chrmName, start, end, minCov, isPrint )
		outAr += regionVals
	# end for i
	# get last bin
	start = nBins * windowSize
	end = chrmLen
	regionVals = processRegion(fastaFile, bamFile, chrmName, start, end, minCov, isPrint)
	outAr += regionVals 
	
	return outAr
	
def processRegion( fastaFile, bamFile, chrm, start, end, minCov, isPrint ):
	
	regionLength = end - start + 1
	#reads = bamFile.fetch(chrm, start, end, multiple_iterators=True)
	reads = bamFile.fetch(chrm, start, end)
	minStart = start
	maxEnd = end
	
	forwardReads = []
	reverseReads = []
	
	for read in reads:
		# read covers minCov of region
		if read.get_overlap(start, end) >= (regionLength * minCov):
			rStart = read.reference_start
			rEnd = read.reference_end
			rLen = read.reference_length
			rSeq = read.query_sequence
			isRev = read.is_reverse
			aRead = (rStart, rSeq)
			if isRev:
				reverseReads += [(rStart, rSeq, rLen)]
			else:
				forwardReads += [(rStart, rSeq, rLen)]
			if rStart < minStart:
				minStart = rStart
			if rEnd > maxEnd:
				maxEnd = rEnd
		# end if read overlap
	# end for read
	refSeq = fastaFile.fetch( chrm, minStart, maxEnd )
	
	#if isPrint:
	#	print( ' ~{:s}:{:d}-{:d}  {:d} forward, {:d} reverse'.format(chrm, start, end, len(forwardReads), len(reverseReads) ) )
	# get read weighted methylation
	weiMeth = computeReads( refSeq, forwardReads, reverseReads, minStart )
	
	# group by level and prep for output
	outAr = groupMethylLevels( weiMeth, start, end )
	return outAr

def computeReads( refSeq, forwardReads, reverseReads, start ):
	outVals = []
	# iterate through forward reads
	for fRead in forwardReads:
		fStart = fRead[0]
		fSeq = fRead[1]
		fLen = fRead[2]
		offset = fStart - start
		fRef = refSeq[ offset: offset+fLen ]
		fZip = zip(fRef, fSeq)
		fVal = computeForwardReadMeth( fZip )
		outVals += [fVal]
	
	# iterate through reverse reads
	for rRead in reverseReads:
		rStart = rRead[0]
		rSeq = rRead[1]
		rLen = rRead[2]
		offset = rStart - start
		rRef = refSeq[ offset: offset+rLen ]
		rZip = zip(rRef, rSeq)
		rVal = computeReverseReadMeth( rZip )
		outVals += [rVal]
	
	return outVals

def computeForwardReadMeth( read ):
	readL = list(read)
	# keep ref 'C' positions
	cPos = list(filter(lambda r: r[0] == 'C', readL))
	# methylated positions -> have 'C'
	mPos = list(filter(lambda r: r[1] == 'C', cPos))
	#print(len(cPos), len(mPos), -1 if len(cPos)== 0 else float(len(mPos)) / float(len(cPos)))
	if len(cPos) == 0:
		return 0
	return float(len(mPos)) / float(len(cPos))
	
def computeReverseReadMeth( read ):
	readL = list(read)
	# keep ref 'G' positions
	cPos = list(filter(lambda r: r[0] == 'G', readL))
	# methylated positions -> have 'C'
	mPos = list(filter(lambda r: r[1] == 'G', cPos))
	#print(len(cPos), len(mPos), -1 if len(cPos)==0 else float(len(mPos)) / float(len(cPos)))
	if len(cPos) == 0:
		return 0
	return float(len(mPos)) / float(len(cPos))

def groupMethylLevels( methylAr, start, end ):
	'''
		output array of tuples with counts per methylation grouping
		currently 10 groups
	'''
	countAr = [0] * 10
	for m in methylAr:
		# m * 100 -> to percent // 10 -> to bin
		mx = m * 100
		i = int( mx // 10 )
		if i == 10:
			i -= 1 # adjust when 100% methylation 
		try:
			countAr[i] += 1
		except IndexError:
			print('ERROR out of range: {:d}, methyl {:g}'.format(i, m))
	# end for m
	
	outAr = []
	total = len(methylAr)
	for i in range(len(countAr)):
		tmp = str(countAr[i]) if total > 0 else 'NA'
		tmp2 = str(float(countAr[i]) / total) if total > 0 else 'NA'
		outAr += [(start, end, i*10, tmp, tmp2)]
	# end for i
	return outAr
	
def writeOutput( outFileStr, outMat, genomeAr, info ):
	outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	headerAr = ['chrm', 'start', 'end', 'methLevel', 'count', 'freq']
	outFile.write( '\t'.join(headerAr) + '\n' )
	
	# loop through chrm
	for i in range(len(genomeAr)):
		chrmName = genomeAr[i][0]
		# loop thru chrm results
		chrmRes = outMat[i]
		for j in range(len(chrmRes)):
			res = chrmRes[j]
			# (0) start, (1) end, (2) mLevel, (3) count [str]
			outStr = '{:s}\t{:d}\t{:d}\t{:d}%\t{:s}\t{:s}\n'.format( chrmName, res[0], res[1], res[2], res[3], res[4] )
			outFile.write( outStr )
		# end for j
	# end for i
	outFile.close()

def parseInputs( argv ):

	outID = None
	isPrint = True
	windowSize = WINDSIZE
	minCov = MINCOV
	startInd = 0
	
	for i in range(len(argv)):
		if argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith('-o='):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-w=' ):
			try:
				windowSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: window size must be integer' )
				exit()
		elif argv[i].startswith('-v='):
			try:
				minCov = float(argv[i][3:])
				startInd += 1
				if minCov > 1.0:
					print( 'ERROR: min cov must be less than 1' )
					exit()
			except ValueError:
				print('ERROR: min cov must be a float')
				exit()
		# end elif
	# end for
	
	fastaFile = argv[startInd]
	bamFile = argv[startInd +1]
	
	processInputs( fastaFile, bamFile, outID, windowSize, minCov, isPrint )	

def printHelp():
	print('Usage: ')


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
