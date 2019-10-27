import sys, math, os, bisect, random, gzip
import pysam
from DmrPos import *
from bisect_func import *

# Usage: region_haplo_allele_v4.py [-q] [-h] [-z] [-rm-snps] [-d] [-b] [-u=unmethyl_val] 
# [-p=min_pos_read][r=min_read]
# [-m=methyl_val] [-o=out_id] [-v=thresh] <pos_file> <bed_file>  <bam_file>

#TRANSITIONS = ['uu', 'uM', 'Mu', 'MM','u.','M.','.u','.M']
#TRANSITIONS = ['uu', 'uM', 'Mu', 'MM']
MTHRESH = 0.35
UNMETHYL = 0.25
METHYL = 0.75
MINREADS=10
MINPOS=2

def processInputs( bamFileStr, dmrPosFileStr, outID, mThresh, debugFile, keepSNPs, unmethylVal, methylVal, minReads, minPos, isBinary, isCompress, isPrint ):
	
	if isPrint:
		print( 'DMR-Position file:', os.path.basename(dmrPosFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'Methylation threshold: {:.3f}'.format(mThresh) )
		print( 'Keep reads with SNPs:', keepSNPs )
		print( 'Max unmethyl value: {:.3f}'.format(unmethylVal) )
		print( 'Min methyl value: {:.3f}'.format(methylVal) )
		print( 'Min reads: {:d}'.format(minReads) )
		print( 'Min positions per read: {:d}'.format(minPos) )
		print( 'Keep only binary reads:', str(isBinary) )
	
	# check BAM file/index
	checkBAM( bamFileStr )

	# determine keep positions in each DMR
	dmrPosAr = readDmrPosFile( dmrPosFileStr )
	if isPrint:
		print( len(dmrPosAr), 'regions found' )
	
	# create pysam bam and fasta objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	
	# do analysis
	outMatrix = []
	debugMatrix = []
	genoMatrix = []
	for dmr in dmrPosAr:
		resDict, debStr = processRegion( bamFile, dmr, mThresh, keepSNPs, minPos, isBinary, isPrint )
		#outMatrix += [resDict]
		debugMatrix += [debStr]
		# format for normal output
		outStr = formatForOutput( dmr, resDict )
		outMatrix += [outStr]
		# call epigenotype and format for output
		genoStr = assignEpigeno( dmr, resDict, unmethylVal, methylVal, minReads )
		genoMatrix += [genoStr]
	
	# close bam and fasta objects
	bamFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
		
	outFileStr = outID + '_dmr_reads-v4.tsv'
	outFileStr += '.gz' if isCompress else ''
	genoFileStr = outFileStr.replace('.tsv', '_geno.tsv')
	if isPrint:
		print( 'Writing output to', outFileStr)
		print( 'Writing genotypes to', genoFileStr )
	
	info = '#from_script: region_haplo_allele_v4.py; dmr_file:{:s}; mThresh:{:.3f}; keep_snps:{:s}; unmethyl_val:{:.3f}; methyl_val:{:.3f}; min_reads:{:d}; min_pos_read:{:d}; is_binary:{:s}\n'.format(os.path.basename(dmrPosFileStr), mThresh, str(keepSNPs), unmethylVal, methylVal, minReads, minPos, str(isBinary))
	outHeader = ['chrm', 'start', 'end', 'label', 'mType', 'count']
	writeOutput( outFileStr, outMatrix, info, outHeader, isCompress )
	
	genoHeader = ['chrm', 'start', 'end', 'label', 'unmethylReads', 'methylReads', 'unmethylLevel', 'methylLevel', 'epigenotype']
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

def readDmrPosFile( dmrPosFileStr ):
	if dmrPosFileStr.endswith('.gz'):
		inFile = gzip.open( dmrPosFileStr, 'rt' )
	else:
		inFile = open( dmrPosFileStr, 'r' )
	
	outAr = []
	
	for line in inFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or lineAr[1].isdigit() == False:
			continue
		chrm = lineAr[0]
		start = int(lineAr[1])
		end = int(lineAr[2])
		label = lineAr[3]
		newRegion = DmrPos( chrm, start, end, label=label )
		posAr = [int(x) for x in lineAr[4].split(',')]
		negAr = [int(x) for x in lineAr[5].split(',')]
		newRegion.setPos( posAr, negAr )
		outAr += [ newRegion ]
	# end for line
	inFile.close()
	return outAr
			
def processRegion( bamFile, dmrRegion, mThresh, keepSNPs, minPos, isBinary, isPrint ):
	#chrm, start, end, label, posPlus, posMinus = dmrRegion
	chrm, start, end, label = dmrRegion.getInfo()
	posPlus = dmrRegion.posAr
	posMinus = dmrRegion.negAr
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
			w = computeReverse( posMinus, rPos, rSeq, keepSNPs, minPos, isBinary )
			ww = encodeReverseRead( posMinus, rPos, rSeq )
			w1 = '- |' + ('' if ww == None else ww) + '|'
		else:
			w = computeForward( posPlus, rPos, rSeq, keepSNPs, minPos, isBinary )
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

def computeForward( posList, readPos, readSeq, keepSNPs, minPos, isBinary ):
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
	if mC + uC < minPos:
		return -1 # not enough positions
	elif isBinary and mC != 0 and uC != 0:
		return -1 # mixed methylation
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

def computeReverse( posList, readPos, readSeq, keepSNPs, minPos, isBinary ):
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
	if mC + uC < minPos:
		return -1 # not enough positions
	elif isBinary and mC != 0 and uC != 0:
		return -1 # mixed methylation
	return float(mC)/ float(mC + uC)

def formatForOutput( dmr, inDict ):
	#chrm, start, end, label, posPlus, posMinus = dmr
	chrm, start, end, label = dmr.getInfo()
	outStr = ''
	
	for x in ['u', 'M', 'x']:
		val = inDict[x]
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:d}\n'.format(chrm, start, end, label, x, val)
	
	return outStr
	
def assignEpigeno( dmr, inDict, minVal, maxVal, minReads ):
	#chrm, start, end, label, posPlus, posMinus = dmr
	chrm, start, end, label = dmr.getInfo()
	
	m = inDict['M']
	u = inDict['u']
	n = m + u
	uLevel = -1 if n < minReads else float(u) / n
	rLevel = -1 if n < minReads else float(m) / n
	if n < minReads:
		#uLevel = -1
		#rLevel = -1
		geno = 'NA'
	elif rLevel <= minVal:
		geno = 'uu' # unemthylated
	elif rLevel <= maxVal:
		geno = 'Mu' # heterozygous methylated
	else:
		geno = 'MM' # methylated
	
	outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:.4f}\t{:.4f}\t{:s}\n'.format(chrm, start, end, label, u, m, uLevel, rLevel, geno )
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
	isBinary = False
	minReads = MINREADS
	minPos = MINPOS
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
		elif argv[i] == '-b':
			isBinary = True
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
		elif argv[i].startswith('-r='):
			try:
				minReads = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: min read parameter invalid...using default', MINREADS)
				minReads = MINREADS
			startInd += 1
		elif argv[i].startswith('-p='):
			try:
				minPos = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: min positions parameter invalid...using default', MINPOS)
				minPos = MINPOS
			startInd += 1
		elif argv[i].startswith('-'):
			print('ERROR: "{:s}" is not a valid parameter\nUse -h to see valid parameters'.format(argv[i]))
			startInd += 1			
		# end elif
	# end for
	
	# check methylVal > unmethylVal
	if unmethylVal > methylVal:
		print('ERROR: Methylated genotype value ({:.3f}) must be greater than unmethylated value ({:.3f})...Use the defaults or check parameters'.format(methylVal, unmethylVal))
		exit()
	
	dmrPosFile = argv[startInd]
	bamFile = argv[startInd +1]
	
	processInputs( bamFile, dmrPosFile, outID, mThresh, debugFile, keepSNPs, unmethylVal, methylVal, minReads, minPos, isBinary, isCompress, isPrint  )
	
def printHelp():
	print('Usage: region_haplo_allele_v3.py [-q] [-h] [-z] [-b] [-rm-snps] [-d] [-u=unmethyl_val] [-m=methyl_val] [-o=out_id] [-v=thresh] <dmr_pos_file>  <bam_file>')


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
