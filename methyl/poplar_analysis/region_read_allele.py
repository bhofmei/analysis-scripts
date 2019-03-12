import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam

# Usage: region_read_allele.py [-q] [-h] [-z] [-o=out_id] <pos_file> <bed_file> <fasta_file> <bam_file> 

def processInputs( fastaFileStr, bamFileStr, bedFileStr, posFileStr, outID, isCompress, isPrint ):
	
	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
	
	# check BAM file/index
	checkBAM( bamFileStr )
	# check FASTA file/index
	checkFASTA( fastaFileStr )
	
	# get DMRs
	dmrAr = readBED( bedFileStr )
	if isPrint:
		print( len(dmrAr), 'regions found' )
	
	# get keep positions - 0-based
	posDict = readPosFile( posFileStr )
	
	# determine keep positions in each DMR
	dmrPosAr = intersectRegionsPos( dmrAr, posDict )
	if isPrint:
		print(len(dmrPosAr), 'regions remaining')
	
	# create pysam bam and fasta objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	fastaFile = pysam.Fastafile( fastaFileStr )
	
	# do analysis
	outMatrix = []
	for dmr in dmrPosAr:
		dmrRes = processRegion( fastaFile, bamFile, dmr, isPrint )
		outMatrix += [dmrRes]
	
	# close bam and fasta objects
	bamFile.close()
	fastaFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_dmr_allele.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: region_read_allele.py; bed_file:{:s}; starting_regions:{:d}; remaining_regions:{:d}\n'.format(os.path.basename(bedFileStr), len(dmrAr), len(dmrPosAr))
	
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
	'''
		Pos file is currently 1-indexed and we need to convert it to 0-indexed
	'''
	
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
		strand = lineAr[2]
		
		if outDict.get(chrm) == None:
			outDict[chrm] = []
		#outDict[chrm] += [pos]
		bisect.insort( outDict[chrm], (pos,strand) )
	# end for line
	posFile.close()
	return outDict

def intersectRegionsPos( dmrAr, posDict ):
	outAr = []
	
	# loop through regions
	for region in dmrAr:
		chrm, start, end, label = region
		
		# pos array for chrm
		if posDict.get(chrm) == None:
			print('WARNING: no valid positions for chrm', chrm, 'but DMR given')
			continue
		posAr = posDict.get(chrm)
		
		# index of first position >= start
		sIndex = bisect_ge( posAr, (start,'') )
		# index of last position < end
		eIndex = bisect_lt( posAr, (end,'') )
		
		# either is None, or eindex <= sIndex , only 1 valid position, so move on
		if eIndex == None or sIndex == None or eIndex <= sIndex:
			continue
		
		# otherwise, get the positions
		cPos = posAr[ sIndex:eIndex+1 ]
		
		# sort the positions into strands
		posPlus = []
		posMinus = []
		for t in cPos:
			if t[1] == '+':
				posPlus += [t[0]]
			elif t[1] == '-':
				posMinus += [t[0]]
		# if only 1 position for both strands, don't keep it
		if len(posPlus) <= 1 and len(posMinus) <= 1:
			continue
		outAr += [(chrm, start, end, label, posPlus, posMinus)]

	# end for region
	return outAr

def bisect_ge(a, x):
	# leftmost (first) INDEX greater than or equal to x
	# return None if not found
	i = bisect.bisect_left(a, x)
	if i != len(a):
		return i
	return None

def bisect_lt(a, x):
	# rightmost (last) INDEX less than x
	# return None if not found
	i = bisect.bisect_left(a, x)
	if i:
		return i-1
	return None
			
def processRegion( fastaFile, bamFile, dmrRegion, isPrint ):
	chrm, start, end, label, posPlus, posMinus = dmrRegion
	# get pairs for comparison
	pairsAr = buildPairs( posPlus, posMinus )
	outStr = ''
	
	if isPrint:
		print(label, '--', len(pairsAr), 'pairs of bases to analyze', pairsAr)
	
	
	
	# loop through pairs
	for i in range(len(pairsAr)):
		# pair is rStart, rEnd (inclusive), strand
		rStart, rEnd, rStrand = pairsAr[i]
		
		# end is exclusive so add 1 to rEnd
		counts = analyzePair( bamFile, chrm, rStrand, rStart, rEnd+1 )
		vals = [ str(x) for x in counts ]
		
		# chrm, start, end, label, pos-pair, uu, uM, Mu, MM, N
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:s}\n'.format(chrm, start, end, label, rStart, rEnd, rStrand, i, '\t'.join(vals))
	# end for i
	return outStr
		
def buildPairs( plusPos, minusPos ):
	outAr = []
	
	for i in range(len(plusPos)-1):
		rStart = plusPos[i]
		rEnd = plusPos[i+1]
		outAr += [(rStart, rEnd, '+')]
		
	for j in range(len(minusPos)-1):
		rStart = minusPos[j]
		rEnd = minusPos[j+1]
		outAr += [(rStart, rEnd, '-')]
	
	outAr.sort() # sort
	return outAr
		

def analyzePair( bamFile, chrm, strand, posA, posB ):
	# posB is exclusive
	
	# get the reads
	reads = bamFile.fetch( chrm, posA, posB )
	
	counts = [0]*6 # uu, uM, Mu, MM, N, R
	typesAr = ['uu', 'uM', 'Mu', 'MM']
	
	for read in reads:

		regionOverlap = read.get_overlap(posA, posB)
		rEnd = read.reference_end
		# must overlap both, so regionOverlap >= posB-posA
		if regionOverlap < (posB-posA):
			continue
		# read doesn't actually overlap posB
		#if posB >= rEnd:
		#	continue
		counts[5] += 1
		rStart = read.reference_start
		rSeq = read.query_sequence
		isRev = read.is_reverse
		
		aType = None
		
		# identify allele type
		if isRev and strand == '-':
			counts[4] += 1
			aType = determineReverseType( rStart, rSeq, posA, posB )
		elif not isRev and strand == '+':
			counts[4] += 1
			aType = determineForwardType( rStart, rSeq, posA, posB )
		
		if aType != None:
			j = typesAr.index(aType)
			counts[j] += 1
	# end for read
	
	return counts

def determineForwardType( fStart, fSeq, posA, posB ):
	nucToType = { 'A':None, 'C':'M','G':None, 'N':None, 'T':'u' }
	offset =  posA - fStart
	endset = (posB-1) - fStart
	try:
		nucA = fSeq[offset]
		nucB = fSeq[endset]
		tA = nucToType[nucA]
		tB = nucToType[nucB]
		#print('+', fStart, posA, posB, offset, endset, posB-posA, len(fSeq))
		#print(fSeq)
		#print(nucA, tA, nucB, tB)
		if tA != None and tB != None:
			return tA+tB
		else:
			return None
	except IndexError:
		# this happens with deletions
		#print(fStart, posA, posB, offset, endset, len(fSeq), fSeq)
		return None

def determineReverseType( rStart, rSeq, posA, posB ):
	nucToType = { 'A':'u', 'C':None,'G':'M', 'N':None, 'T':None }
	offset =  posA - rStart
	endset = (posB-1) - rStart
	try:
		nucA = rSeq[offset]
		nucB = rSeq[endset]
		tA = nucToType[nucA]
		tB = nucToType[nucB]
		#print('-', rStart, posA, posB, offset, endset, posB-posA, len(rSeq))
		#print(rSeq)
		#print(nucA, tA, nucB, tB)
		if tA != None and tB != None:
			return tA+tB
		else:
			return None
	except IndexError:
		# this happens with deletions
		#print(rStart, posA, posB, offset, endset, len(rSeq), rSeq)
		return None

def writeOutput( outFileStr, outMat, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	headerAr = ['chrm', 'start', 'end', 'label', 'posA', 'posB','strand', 'i','uu', 'uM', 'Mu', 'MM', 'N','R']
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
		elif argv[i].startswith('-'):
			print('ERROR: "{:s}" is not a valid parameter\nUse -h to see valid parameters'.format(argv[i]))
			
		# end elif
	# end for
	
	posFile = argv[startInd]
	bedFile = argv[startInd+1]
	fastaFile = argv[startInd+2]
	bamFile = argv[startInd +3]
	
	
	processInputs( fastaFile, bamFile, bedFile, posFile, outID, isCompress, isPrint )
	
def printHelp():
	print('Usage: python3 region_read_allele.py [-q] [-h] [-o=out_id] <pos_file> <bed_file> <fasta_file> <bam_file>')


if __name__ == "__main__":
	if len(sys.argv) < 5 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
