import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam

MINCOV=0.5
METHTYPE='C'

# Usage: python3 wei_meth_for_region_reads3.py [-h] [-q] [-z] [-o=outID] [-v=min_overlap] <pos_file> <fasta_file> <bam_file> <bed_file>

def processInputs( fastaFileStr, bamFileStr, bedFileStr, posFileStr, outID, minCov, isCompress, isPrint ):
	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'BED file:', bedFileStr )
		print( 'Min percent overlap:', minCov )
	
	# check BAM file/index
	checkBAM( bamFileStr )
	# check FASTA file/index
	checkFASTA( fastaFileStr )
	
	# get DMR regions
	dmrAr = readBED( bedFileStr )
	
	# get keep positions - 0-based
	posDict = readPosFile( posFileStr )
	
	# create pysam file objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	fastaFile = pysam.Fastafile( fastaFileStr )
	
	# process
	if isPrint:
		print('Analyzing {:d} DMRs'.format( len(dmrAr)))
		
	outMatrix = []
	for dmr in dmrAr:
		dmrRes = processRegion( fastaFile, bamFile, posDict, dmr, minCov, isPrint )
		outMatrix += [dmrRes]
	
	bamFile.close()
	fastaFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_dmr_mReads3.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: wei_meth_for_region_reads3.py; bed_file:{:s}; min_cov:{:g}\n'.format(os.path.basename(bedFileStr), minCov)
	
	writeOutput( outFileStr, outMatrix, info, isCompress )
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

def readBED( bedFileStr ):
	
	bedFile = open(bedFileStr, 'r' )
	outAr = []
	count = 1
	
	for line in bedFile:
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) start (2) end (3) name?
		if len(lineAr) < 3:
			continue
		chrm = lineAr[0]
		start = int(lineAr[1])
		end = int(lineAr[2])
		dmrName = ( 'DMR-'+count if len(lineAr) < 4 else lineAr[3] )
		outAr += [(chrm, start, end, dmrName)]
		count += 1
	# end for line
	bedFile.close()
	return outAr

def readPosFile( posFileStr ):
	
	if posFileStr.endswith('.gz'):
		posFile = gzip.open( posFileStr, 'rt' )
	else:
		posFile = open( posFileStr, 'r' )
	outDict = {}
	
	for line in posFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or lineAr[1].isdigit() == False:
			continue
		chrm = lineAr[0]
		pos = int(lineAr[1]) - 1
		
		if outDict.get(chrm) == None:
			outDict[chrm] = []
		outDict[chrm] += [pos]
	# end for line
	posFile.close()
	return outDict
	
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
		
def processRegion( fastaFile, bamFile, posDict, dmrRegion, minCov, isPrint ):
	
	chrm, start, end, label = dmrRegion
	posAr = []
	if posDict.get(chrm) == None:
		print('WARNING: no valid positions for chrm', chrm, 'but DMR given')
	else:
		posAr = posDict[chrm]
	regionLength = end - start + 1
	reads = bamFile.fetch(chrm, start, end)
	
	forwardReads = []
	reverseReads = []
	
	for read in reads:
		# read covers minCov of region
		# if readLen < regionLength, keep reads with 50% overlap
		# if readLen > regionLength, keep 
		regionOverlap = read.get_overlap(start, end)
		readLength = read.query_length
		'''print( '***** readLength, regionLength, regionOverlap', readLength, regionLength, regionOverlap)
		print('readLength <= regionLength', readLength <= regionLength )
		print('regionOverlap >= (readLength * minCov) ', regionOverlap >= (readLength * minCov) )
		print( 'readLength > regionLength', readLength > regionLength)
		print('regionOverlap >= (regionLength * minCov)', regionOverlap >= (regionLength * minCov))'''
		if (minCov == 0) or ( readLength <= regionLength and regionOverlap >= (readLength * minCov) ) or (readLength > regionLength and regionOverlap >= (regionLength * minCov)):
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
		# end if read overlap
	# end for read
	refSeq = fastaFile.fetch( chrm, start-2, end+2 )
	
	# get the positions we want - if none, mo
	posKeep = offsetPosList( posAr, start, end )
	
	if len(posAr) == 0:
		return  formatMethylLevels( [], dmrRegion )
	
	# get read weighted methylation
	weiMeth = computeReads( refSeq, forwardReads, reverseReads, start, end, posKeep )
	
	# prep for output
	outStr = formatMethylLevels( weiMeth, dmrRegion )
	return outStr

def offsetPosList( posAr, start, end ):
	# start index is left-most greater or equal to start
	i = bisect.bisect_left(posAr, start)
	# end index is right-most less than end (end not inclusive)
	j = bisect.bisect_left(posAr, end)
	if i >= j:
		return []
	#print( '*', start, end, i, j)
	#print('*', posAr[i:j])
	outAr = [ posAr[k] - start for k in range(i, j) ]
	#print('*', outAr)
	return outAr

def bisectContains( a, x ):
	i = bisect.bisect_left( a, x )
	return (i != len(a) and a[i] == x)

def computeReads( refSeq, forwardReads, reverseReads, start, end, posKeep ):
	outVals = []
	regionLen = end - start + 1
	# iterate through forward reads
	for fRead in forwardReads:
		fStart = fRead[0]
		fSeq = fRead[1]
		fLen = fRead[2]
		offset, endset, frontPadding = readPortion( start, end, fStart, fLen )
		fSeqs = 'N' * frontPadding + fSeq[offset:endset]
		fRef = refSeq[ 2: ]
		fVal = computeForwardReadMeth( fRef, fSeqs, posKeep )
		outVals += [fVal]
	
	# iterate through reverse reads
	for rRead in reverseReads:
		rStart = rRead[0]
		rSeq = rRead[1]
		rLen = rRead[2]
		offset, endset, frontPadding = readPortion( start, end, rStart, rLen )
		rSeqs = 'N' * frontPadding + rSeq[offset:endset]
		rRef = refSeq[ : -2 ]
		rVal = computeReverseReadMeth( rRef, rSeqs, posKeep )
		outVals += [rVal]
	
	return outVals
	
def readPortion( start, end, readStart, readLen ):
	if readStart <= start:
		offset = start - readStart
		frontPadding = 0
		if (readStart + readLen) < end:
			# only front part
			endset = readLen
		else:
			# middle
			endset = end - readStart
	else:
		# only back part
		frontPadding = readStart - start
		offset = 0
		endset = end - readStart
	#print(start, end, readStart, readLen, offset, endset, frontPadding)
	return (offset, endset, frontPadding)

def computeForwardReadMeth( sRef, sSeq, posKeep ):

	cgPos = []
	chgPos = []
	chhPos = []
	readLen = len(sSeq)
	# 	print('*****')
	# 	print(sRef)
	# 	print(sSeq)
	# 	print(posKeep)
	
	#for i in range(min(len(sRef)-2, len(sSeq))):
	for i in posKeep:
		# don't include, not keep or not in read
		#if bisectContains( posKeep, i ) == False or sSeq[i] == 'N':
		if i >= readLen:
			break
		if sSeq[i] == 'N':
			continue
		if sRef[i] == 'C' and sRef[i+1] == 'G':
			cgPos += [sSeq[i]]
		elif sRef[i] == 'C' and sRef[i+2] == 'G':
			chgPos += [sSeq[i]]
		elif sRef[i] == 'C':
			chhPos += [sSeq[i]]
	# end for
	mcgPos = list(filter(lambda r: r == 'C', cgPos))
	mchgPos = list(filter(lambda r: r == 'C', chgPos))
	mchhPos = list(filter(lambda r: r == 'C', chhPos))
	
	nCg = len(cgPos)
	nChg = len(chgPos)
	nChh = len(chhPos)
	nC = nCg + nChg + nChh
	
	mCg = len(mcgPos)
	mChg = len(mchgPos)
	mChh = len(mchhPos)
	mC = mCg + mChg + mChh
	
	# compute weights
	wC = -1 if nC == 0 else float(mC) / nC
	wCg = -1 if nCg == 0 else float(mCg) / nCg
	wChg = -1 if nChg == 0 else float(mChg) / nChg
	wChh = -1 if nChh == 0 else float(mChh) / nChh
	#print('--', mCg, nCg)
	return (wC, mC, nC, wCg, mCg, nCg, wChg, mChg, nChg, wChh, mChh, nChh)
	
def computeReverseReadMeth( sRef, sSeq, posKeep ):
	
	cgPos = []
	chgPos = []
	chhPos = []
	readLen = len(sSeq)
	# 	print('*****')
	# 	print(sRef)
	# 	print('  '+sSeq)
	# 	print(posKeep)
	
	#for i in range(min(len(sRef)-2, len(sSeq))):
	for i in posKeep:
		#if bisectContains( posKeep, i ) == False or sSeq[i] == 'N':
		if i >= readLen:
			break
		if sSeq[i] == 'N':
			continue
		if sRef[i+2] == 'G' and sRef[i+1] == 'C':
			cgPos += [sSeq[i]]
		elif sRef[i+2] == 'G' and sRef[i] == 'C':
			chgPos += [sSeq[i]]
		elif sRef[i+2] == 'G':
			chhPos += [sSeq[i]]
	# end for
	mcgPos = list(filter(lambda r: r == 'G', cgPos))
	mchgPos = list(filter(lambda r: r == 'G', chgPos))
	mchhPos = list(filter(lambda r: r == 'G', chhPos))
		
	nCg = len(cgPos)
	nChg = len(chgPos)
	nChh = len(chhPos)
	nC = nCg + nChg + nChh
	
	mCg = len(mcgPos)
	mChg = len(mchgPos)
	mChh = len(mchhPos)
	mC = mCg + mChg + mChh
	
	# compute weights
	wC = -1 if nC == 0 else float(mC) / nC
	wCg = -1 if nCg == 0 else float(mCg) / nCg
	wChg = -1 if nChg == 0 else float(mChg) / nChg
	wChh = -1 if nChh == 0 else float(mChh) / nChh
	#print('--', mCg, nCg)
	return (wC, mC, nC, wCg, mCg, nCg, wChg, mChg, nChg, wChh, mChh, nChh)

def formatMethylLevels( methylAr, dmrRegion ):
	'''
		output string
	'''
	chrm, start, end, label = dmrRegion
	outData = ''
	if len(methylAr) == 0:
		valStr = ['NA' for x in range(12)]
		return '{:s}\t{:d}\t{:d}\t{:s}\tNA\t{:s}\n'.format( chrm, start, end, label,  '\t'.join(valStr) )
		
	for i in range(len(methylAr)):
		# outStr = chrm, start, end, label, i, level, mCs, tCs, mCGs, tCGs, mCHGs, tCHGs, mCHHs, tCHHs
		val = methylAr[i]
		valStr = ['NA' if x == -1 else str(x) for x in val]
		outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:s}\n'.format( chrm, start, end, label, i, '\t'.join(valStr) )
		outData += outStr
	return outData
	
def writeOutput( outFileStr, outMat, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	headerAr = ['chrm', 'start', 'end', 'label', 'i', 'wmC', 'mC','tC', 'wmCG', 'mCG', 'tCG', 'wmCHG', 'mCHG', 'tCHG', 'wmCHH', 'mCHH', 'tCHH']
	outFile.write( '\t'.join(headerAr) + '\n' )
	
	# loop through regions
	for x in outMat:
		outFile.write(x)
	# end for x
	outFile.close()

def parseInputs( argv ):

	outID = None
	isCompress = False
	isPrint = True
	minCov = MINCOV
	startInd = 0
	
	for i in range(len(argv)):
		if argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i].startswith('-o='):
			outID = argv[i][3:]
			startInd += 1
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
	
	posFile = argv[startInd]
	fastaFile = argv[startInd+1]
	bamFile = argv[startInd +2]
	bedFile = argv[startInd+3]
	
	processInputs( fastaFile, bamFile, bedFile, posFile, outID, minCov, isCompress, isPrint )	

def printHelp():
	print('Usage: python3 wei_meth_for_region_reads.py [-h] [-q] [-z] [-o=outID] [-v=min_overlap] [-m=meth_type] <pos_file> <fasta_file> <bam_file> <bed_file>')


if __name__ == "__main__":
	if len(sys.argv) < 5 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
