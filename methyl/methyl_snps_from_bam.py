import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam #, pybedtools

# Usage: python methyl_snps_from_bam.py [-h] [-q] [-z] [-v=min_overlap] [-o=out_id] <fasta_file> <bed_file> <bam_file>

MINOVER=1

def processInputs( fastaFileStr, bamFileStr, bedFileStr, outID, minOverlap, isCompress, isPrint ):
	if isPrint:
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'BED file:', bedFileStr )
		print( 'Min basepair overlap:', minOverlap )
	
	# check BAM file/index
	checkBAM( bamFileStr )
	# check FASTA file/index
	checkFASTA( fastaFileStr )
	
	# get DMR regions
	regionAr = readBED( bedFileStr )
	if isPrint:
		print('{:d} regions'.format(len(regionAr)))
	
	# create pysam file objects
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	fastaFile = pysam.Fastafile( fastaFileStr )
	
	outMatrix = []
	for i in range(len(regionAr)):
		region = regionAr[i]
		if isPrint and (i % 1000 == 0):
			print('Processed {:d} regions...'.format(i))
			
		res = processRegion( fastaFile, bamFile, region, minOverlap )
		outMatrix += [res]
	
	bamFile.close()
	fastaFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_meReads.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: methyl_snps_from_bam.py; bed_file:{:s}; min_overlap:{:g}\n'.format(os.path.basename(bedFileStr), minOverlap)
	
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
		name = ( 'region-'+count if len(lineAr) < 4 else lineAr[3] )
		outAr += [(chrm, start, end, name)]
		count += 1
	# end for line
	bedFile.close()
	return outAr

def processRegion( fastaFile, bamFile, region, minOverlap ):
	
	chrm, start, end, label = region
	reads = bamFile.fetch(chrm, start, end)
	refSeq = fastaFile.fetch( chrm, start, end )
	
	outStr = '## {:s}\t{:s}:{:d}-{:d}\n'.format(label, chrm, start, end)
	outStr += '  ' + refSeq +'\n'
	
	for read in reads:
		regionOverlap = read.get_overlap(start, end)
		if regionOverlap >= minOverlap :
			rPos = read.get_reference_positions(full_length=True)
			rSeq = read.query_sequence
			isRev = read.is_reverse
			if isRev:
				w = encodeReverse( start, refSeq, rPos, rSeq )
			else:
				w = encodeForward( start, refSeq, rPos, rSeq )
			outStr += ('' if w == None else w)
		# end if read overlap
	# end for read
	
	# add an extra blank line
	outStr += '\n'
	return outStr

def encodeForward( rStart, refSeq, readPos, readSeq ):
	outStr = '+ '
	found = 0
	
	for i in range(len(refSeq)):
		j = bisect_index(readPos, rStart + i)
		if j == None:
			outStr += ' '
		else:
			found += 1
			ref = refSeq[i]
			nuc = readSeq[j]
			if ref == 'C' and nuc == 'C':
				outStr += 'M' # methylated
			elif ref == 'C' and nuc == 'T':
				outStr += 'u' # unmethylated
			elif ref != nuc:
				outStr += nuc # snp
			else:
				outStr += '_'
	# end for
	outStr += '\n'
	return (outStr if found > 1 else None)

def encodeReverse( rStart, refSeq, readPos, readSeq ):
	outStr = '- '
	found = 0
	
	for i in range(len(refSeq)):
		j = bisect_index(readPos, rStart + i)
		if j == None:
			outStr += ' '
		else:
			found += 1
			ref = refSeq[i]
			nuc = readSeq[j]
			rComp = baseComplement(nuc)
			if ref == 'G' and nuc == 'G':
				outStr += 'M'
			elif ref == 'G' and nuc == 'A':
				outStr += 'u'
			elif ref != nuc:
				outStr += rComp
			else:
				outStr += '_'
	# end for
	outStr += '\n'
	return (outStr if found > 1 else None)

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
			
def baseComplement( n ):
	baseAr = ['A','C','G','T']
	revAr = ['T','G','C','A']
	try:
		ind = baseAr.index( n.upper() )
		return revAr[ind]
	except ValueError:
		return 'X'

def writeOutput( outFileStr, outMat, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info + '\n')
	
	# loop through regions
	for x in outMat:
		outFile.write(x)
	# end for x
	outFile.close()

def parseInputs( argv ):
	outID = None
	minOverlap = MINOVER
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
		elif argv[i].startswith('-v='):
			try:
				minOverlap = int(argv[i][3:])
			except ValueError:
				print('WARNING: min overlap parameter invalid...using default', MINOVER)
				minOverlap = MINOVER
			startInd += 1
		elif argv[i].startswith('-'):
			print('ERROR: "{:s}" is not a valid parameter\nUse -h to see valid parameters'.format(argv[i]))
	# end for
	
	fastaFileStr = argv[startInd]
	bedFileStr = argv[startInd+1]
	bamFileStr = argv[startInd+2]
	
	processInputs( fastaFileStr, bamFileStr, bedFileStr, outID, minOverlap, isCompress, isPrint )	

def printHelp():
	print('Usage: python methyl_snps_from_bam.py [-h] [-q] [-z] [-v=min_overlap] [-o=out_id] <fasta_file> <bed_file> <bam_file>')
	print()
	print('Required:')
	print('fasta_file\treference sequence used for mapping')
	print('bed_file\tlist of regions to get reads for')
	print('bam_file\tmapped methylC-seq reads')
	print('Optional:')
	print('-h\t\tprint this help message')
	print('-q\t\tquiet; do not print progress')
	print('-z\t\tcompress output file with gzip')
	print('-o=out_id\toutput file identifier; default is bam file name')
	print('-v=min_overlap\tmin number of bp a read must overlap a region to be included [default 1]')
	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
