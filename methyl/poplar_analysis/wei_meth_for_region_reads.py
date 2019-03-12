import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam

MINCOV=0.5
METHTYPE='C'

# Usage: python3 wei_meth_for_region_reads.py [-h] [-q] [-z] [-o=outID] [-v=min_overlap] <fasta_file> <bam_file> <bed_file>

def processInputs( fastaFileStr, bamFileStr, bedFileStr, outID, minCov, isCompress, isPrint ):
	if isPrint:
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
	
	# create pysam file objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	fastaFile = pysam.Fastafile( fastaFileStr )
	
	# process
	if isPrint:
		print('Analyzing {:d} DMRs'.format( len(dmrAr)))
		
	outMatrix = []
	for dmr in dmrAr:
		dmrRes = processRegion( fastaFile, bamFile, dmr, minCov, isPrint )
		outMatrix += [dmrRes]
	
	bamFile.close()
	fastaFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_dmr_mReads.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: wei_meth_for_region_reads.py; bed_file:{:s}; min_cov:{:g}\n'.format(os.path.basename(bedFileStr), minCov)
	
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
	
def processRegion( fastaFile, bamFile, dmrRegion, minCov, isPrint ):
	
	chrm, start, end, label = dmrRegion
	regionLength = end - start + 1
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
	refSeq = fastaFile.fetch( chrm, minStart-2, maxEnd+2 )
	
	# get read weighted methylation
	weiMeth = computeReads( refSeq, forwardReads, reverseReads, minStart-2 )
	
	# prep for output
	outStr = formatMethylLevels( weiMeth, dmrRegion )
	return outStr

def computeReads( refSeq, forwardReads, reverseReads, start ):
	outVals = []
	# iterate through forward reads
	for fRead in forwardReads:
		fStart = fRead[0]
		fSeq = fRead[1]
		fLen = fRead[2]
		offset = fStart - start+2
		fRef = refSeq[ offset: offset+fLen+2 ]
		fVal = computeForwardReadMeth( fRef, fSeq )
		outVals += [fVal] if fVal != -1 else []
	
	# iterate through reverse reads
	for rRead in reverseReads:
		rStart = rRead[0]
		rSeq = rRead[1]
		rLen = rRead[2]
		offset = rStart - start
		rRef = refSeq[ offset: offset+rLen+2 ]
		rVal = computeReverseReadMeth( rRef, rSeq )
		outVals += [rVal] if rVal != -1 else []
	
	return outVals

def computeForwardReadMeth( sRef, sSeq ):

	cgPos = []
	chgPos = []
	chhPos = []
	
	for i in range(min(len(sRef)-2, len(sSeq))):
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
	
	nC = len(cgPos) + len(chgPos) + len(chhPos)
	if nC == 0:
		return -1
	mC = len(mcgPos) + len(mchgPos) + len(mchhPos)
	k = float(mC) / float(nC)
	return (k, mC, nC, len(mcgPos), len(cgPos), len(mchgPos), len(chgPos), len(mchhPos), len(chhPos))
	
def computeReverseReadMeth( sRef, sSeq ):
	
	cgPos = []
	chgPos = []
	chhPos = []
	for i in range(min(len(sRef)-2, len(sSeq))):
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
		
	nC = len(cgPos) + len(chgPos) + len(chhPos)
	if nC == 0:
		return -1
	mC = len(mcgPos) + len(mchgPos) + len(mchhPos)
	k = float(mC) / float(nC)
	return (k, mC, nC, len(mcgPos), len(cgPos), len(mchgPos), len(chgPos), len(mchhPos), len(chhPos))

def formatMethylLevels( methylAr, dmrRegion ):
	'''
		output string
	'''
	chrm, start, end, label = dmrRegion
	outData = ''
	for i in range(len(methylAr)):
		# outStr = chrm, start, end, label, i, level, mCs, tCs, mCGs, tCGs, mCHGs, tCHGs, mCHHs, tCHHs
		val = methylAr[i]
		valStr = [str(x) for x in val[1:]]
		outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:f}\t{:s}\n'.format( chrm, start, end, label, i, val[0], '\t'.join(valStr) )
		outData += outStr
	return outData
	
def writeOutput( outFileStr, outMat, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	headerAr = ['chrm', 'start', 'end', 'label', 'i', 'methLevel', 'mCs','tCs', 'mCGs', 'tCGs', 'mCHGs', 'tCHGs', 'mCHHs', 'tCHHs']
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
	methType = METHTYPE
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
		elif argv[i].startswith('-m='):
			mType = argv[i][3:].upper()
			if mType not in ['C', 'CG', 'CHG', 'CHH', 'CNN']:
				print('WARNING: invalid methylation type...using default', METHTYPE)
			else:
				methType = 'C' if mType == 'CNN' else mType
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
	
	fastaFile = argv[startInd]
	bamFile = argv[startInd +1]
	bedFile = argv[startInd+2]
	
	processInputs( fastaFile, bamFile, bedFile, outID, minCov, isCompress, isPrint )	

def printHelp():
	print('Usage: ')


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
