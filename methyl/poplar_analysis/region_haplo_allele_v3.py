import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam

# Usage: region_haplo_allele_v3.py [-q] [-h] [-z] [-rm-snps] [-d] [-u=unmethyl_val] [-m=methyl_val] [-o=out_id] [-v=thresh] <pos_file> <bed_file>  <bam_file>

#TRANSITIONS = ['uu', 'uM', 'Mu', 'MM','u.','M.','.u','.M']
#TRANSITIONS = ['uu', 'uM', 'Mu', 'MM']
MTHRESH = 0.35
UNMETHYL = 0.25
METHYL = 0.75

def processInputs( bamFileStr, bedFileStr, posFileStr, outID, mThresh, debugFile, keepSNPs, unmethylVal, methylVal, isCompress, isPrint ):
	
	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
		print( 'Methylation threshold: {:.3f}'.format(mThresh) )
		print( 'Keep reads with SNPs:', keepSNPs )
		print( 'Max unmethyl value: {:.3f}'.format(unmethylVal) )
		print( 'Min methyl value: {:.3f}'.format(methylVal) )
	
	# check BAM file/index
	checkBAM( bamFileStr )
	
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
	
	# do analysis
	outMatrix = []
	debugMatrix = []
	genoMatrix = []
	for dmr in dmrPosAr:
		resDict, debStr = processRegion( bamFile, dmr, mThresh, keepSNPs, isPrint )
		#outMatrix += [resDict]
		debugMatrix += [debStr]
		# format for normal output
		outStr = formatForOutput( dmr, resDict )
		outMatrix += [outStr]
		# call epigenotype and format for output
		genoStr = assignEpigeno( dmr, resDict, unmethylVal, methylVal )
		genoMatrix += [genoStr]
	
	# close bam and fasta objects
	bamFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
		
	outFileStr = outID + '_dmr_reads-v3.tsv'
	outFileStr += '.gz' if isCompress else ''
	genoFileStr = outFileStr.replace('.tsv', '_geno.tsv')
	if isPrint:
		print( 'Writing output to', outFileStr)
		print( 'Writing genotypes to', genoFileStr )
	
	info = '#from_script: region_haplo_allele_v3.py; bed_file:{:s}; starting_regions:{:d}; remaining_regions:{:d}; mThresh:{:.3f}; keep_snps:{:s}; unmethyl_val:{:.3f}; methyl_val:{:.3f}\n'.format(os.path.basename(bedFileStr), len(dmrAr), len(dmrPosAr), mThresh, str(keepSNPs), unmethylVal, methylVal)
	outHeader = ['chrm', 'start', 'end', 'label', 'mType', 'count']
	writeOutput( outFileStr, outMatrix, info, outHeader, isCompress )
	
	genoHeader = ['chrm', 'start', 'end', 'label', 'unmethylReads', 'methylReads', 'epigenotype']
	writeOutput( genoFileStr, genoMatrix, info, genoHeader, isCompress )
	
	if debugFile:
		#debugFileStr = outID + '_dmr_reads-v3_debug.tsv'
		#debugFileStr += '.gz' if isCompress else ''
		debugFileStr = outFileStr.replace('.tsv', '_debug.tsv')
		debugHeader = None
		if isPrint:
			print( 'Writing debug output to', debugFileStr )
		writeOutput( debugFileStr, debugMatrix, info, debugHeader, isCompress )
	
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

def bisect_index(a, x):
	try:
		i = bisect.bisect_left(a, x)
		if i != len(a) and a[i] == x:
			return i
		return None
	except TypeError:
		if x in a:
			return a.index(x)
		else:
			return None
			
def processRegion( bamFile, dmrRegion, mThresh, keepSNPs, isPrint ):
	chrm, start, end, label, posPlus, posMinus = dmrRegion
	outDict = {'u': 0, 'M': 0, 'x': 0}
	
	# get reads
	reads = bamFile.fetch(chrm, start, end)
	outDebug = '## {:s}\n'.format( label)
	# loop through reads
	for read in reads:
		rPos = read.get_reference_positions(full_length=True)
		rSeq = read.query_sequence
		isRev = read.is_reverse
		
		# determine if read is 'methylated' or 'unmethylated'
		if isRev:
			w = computeReverse( posMinus, rPos, rSeq, keepSNPs )
			ww = encodeReverseRead( posMinus, rPos, rSeq)
			w1 = '- |' + ('' if ww == None else ww) + '|'
		else:
			w = computeForward( posPlus, rPos, rSeq, keepSNPs )
			ww = encodeForwardRead( posPlus, rPos, rSeq) 
			w1 = '+ |' + ('' if ww == None else ww) + '|'
		
		# add to dict
		if w == -1:
			w2='x'
			outDict['x'] += 1
		elif w < mThresh:
			w2 = 'u'
			outDict['u'] += 1
		else:
			w2 = 'M'
			outDict['M'] += 1
		outDebug += '{:s} {:.3f} {:s}\n'.format(w1, w, w2)
	# end for read
	
	#return outStr, outDebug
	return outDict, outDebug

def encodeForwardRead( posList, readPos, readSeq ):
	outStr = ''
	found = 0
	for p in posList:
		j = bisect_index(readPos, p)
		if j == None:
			outStr += ' '
		else:
			found += 1
			nuc = readSeq[j]
			if nuc  == 'C':
				outStr += 'M'
			elif nuc == 'T':
				outStr += 'u'
			else:
				outStr += '.'
	# end for p
	return (outStr if found > 1 else None)

def computeForward( posList, readPos, readSeq, keepSNPs ):
	uC = 0
	mC = 0
	#print('+:', posList)
	for p in posList:
		j = bisect_index(readPos, p)
		if j == None:
			continue
		nuc = readSeq[j]
		#print('+', p, j, nuc)
		if nuc == 'C':
			mC += 1 # meth
		elif nuc == 'T':
			uC += 1 # unmeth
		elif not keepSNPs:
			return -1 # SNP and we are throwing away those reads
	# end for p
	if mC + uC < 1:
		return -1 # not enough positions
	return float(mC)/ float(mC + uC)
	
def encodeReverseRead( posList, readPos, readSeq ):
	outStr = ''
	found = 0
	for p in posList:
		j = bisect_index(readPos, p)
		if j == None:
			outStr += ' '
		else:
			found += 1
			nuc = readSeq[j]
			if nuc  == 'G':
				outStr += 'M'
			elif nuc == 'A':
				outStr += 'u'
			else:
				outStr += '.'
	# end for p
	return (outStr if found > 1 else None)

def computeReverse( posList, readPos, readSeq, keepSNPs ):
	uC = 0
	mC = 0
	#print('-:', posList)
	for p in posList:
		j = bisect_index(readPos, p)
		if j == None:
			continue
		nuc = readSeq[j]
		#print( '-', p, j, nuc )
		if nuc == 'G':
			mC += 1 # meth
		elif nuc == 'A':
			uC += 1 # unmeth
		elif not keepSNPs:
			return -1 # SNP and we are throwing away those reads
	# end for p
	#print('--', mC, uC)
	if mC + uC < 1:
		return -1 # not enough positions
	return float(mC)/ float(mC + uC)

def formatForOutput( dmr, inDict ):
	chrm, start, end, label, posPlus, posMinus = dmr
	outStr = ''
	
	for x in ['u', 'M', 'x']:
		val = inDict[x]
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:d}\n'.format(chrm, start, end, label, x, val)
	
	return outStr
	
def assignEpigeno( dmr, inDict, minVal, maxVal ):
	chrm, start, end, label, posPlus, posMinus = dmr
	
	m = inDict['M']
	u = inDict['u']
	n = m + u
	uLevel = float(u) / n
	rLevel = float(m) / n
	if rLevel <= minVal:
		geno = 'uu' # unemthylated
	elif rLevel <= maxVal:
		geno = 'Mu' # heterozygous methylated
	else:
		geno = 'MM' # methylated
	
	outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:.4f}\t{:.4f}\t{:s}\n'.format(chrm, start, end, label, uLevel, rLevel, geno )
	return outStr

def writeOutput( outFileStr, outMat, info, headerAr, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	
	if headerAr != None and len(headerAr) > 0:
		outFile.write( '\t'.join(headerAr) + '\n' )
	
	# loop through regions
	for x in outMat:
		outFile.write(x)
	# end for x
	outFile.close()

def parseInputs( argv ):
	outID = None
	mThresh = MTHRESH
	unmethylVal = UNMETHYL
	methylVal = METHYL
	keepSNPs = True
	debugFile = False
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
		elif argv[i] == '-d':
			debugFile = True
			startInd += 1
		elif argv[i] == '-rm-snps':
			keepSNPs = False
			startInd += 1
		elif argv[i].startswith('-u='):
			try:
				unmethylVal = float(argv[i][3:])
				unmethylVal = ( unmethylVal / 100.0 if unmethylVal > 1 else unmethylVal )
			except ValueError:
				print('Unmethylated genotype value parameter is invalid...using default', UNMETHYL)
				unmethylVal = UNMETHYL
			startInd += 1
		elif argv[i].startswith('-m='):
			try:
				methylVal = float(argv[i][3:])
				methylVal = ( unmethylVal / 100.0 if methylVal > 1 else methylVal )
			except ValueError:
				print('Methylated genotype value parameter is invalid...using default', METHYL)
				methylVal = METHYL
			startInd += 1
		elif argv[i].startswith('-o='):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-v='):
			try:
				mThresh = float(argv[i][3:])
				mThresh = (mThresh / 100.0 if mThresh > 1 else mThresh)
			except ValueError:
				print('WARNING: methylation threshold parameter invalid...using default', MTHRESH)
				mThresh = MTHRESH
			startInd += 1
		elif argv[i].startswith('-'):
			print('ERROR: "{:s}" is not a valid parameter\nUse -h to see valid parameters'.format(argv[i]))
			
		# end elif
	# end for
	
	# check methylVal > unmethylVal
	if unmethylVal > methylVal:
		print('ERROR: Methylated genotype value ({:.3f}) must be greater than unmethylated value ({:.3f})...Use the defaults or check parameters'.format(methylVal, unmethylVal))
		exit()
	
	posFile = argv[startInd]
	bedFile = argv[startInd+1]
	bamFile = argv[startInd +2]
	
	processInputs( bamFile, bedFile, posFile, outID, mThresh, debugFile, keepSNPs, unmethylVal, methylVal,  isCompress, isPrint )
	
def printHelp():
	print('Usage: region_haplo_allele_v3.py [-q] [-h] [-z] [-rm-snps] [-d] [-u=unmethyl_val] [-m=methyl_val] [-o=out_id] [-v=thresh] <pos_file> <bed_file>  <bam_file>')


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
