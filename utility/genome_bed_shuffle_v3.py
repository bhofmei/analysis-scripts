import sys, math, glob, multiprocessing, subprocess, os, bisect, random, re
import pysam, pybedtools

# Usage: python genome_bed_shuffle_v3.py [-h] [-q] [-k] [-s=seed_val] 
# [-n=max_tries] [-o=out_id] [-t=pos_file] [-p=pos_threshold] [-b=bam_files]
# [-r=min_read_cov] [-i=include_file] <fasta_file> <bed_file>

NTRIES=1000
POSTHRESH = -1
NREADS=5

def processInputs( bedFileStr, fastaFileStr, posFileStr, bamFileStr, minReads, maxTries, outId, seedVal, posThresh, includeFileStr, isKeep, isPrint ):

	isPos = (posFileStr != None)
	isBam = (bamFileStr != None)
	
	if outId == None:
		tmp = os.path.basename( bedFileStr )
		outId = re.sub( "\.bed", "", tmp, flags=re.IGNORECASE )
	outFileStr = outId + '_equiv' + '.bed'
	
	if isBam:
		bamFileAr = readBamList( bamFileStr )
	
	# check include file
	if includeFileStr != None and os.path.isfile(includeFileStr) == False:
		print( 'WARNING: include file {:s} does not exist...ignoring it'.format(os.path.basename(includeFileStr)))
		includeFileStr = None
		
	if isPrint:
		print( 'Max Tries:', maxTries )
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		if isPos:
			print( 'Position Threshold:', posThresh)
			print( 'Position file:', os.path.basename(posFileStr) )
		print( 'Starting seed:', seedVal )
		if isBam:
			print( 'BAM files:', len(bamFileAr) )
			print( 'Min read coverage:', minReads )
		if includeFileStr != None:
			print( 'Include regions file:', os.path.basename(includeFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
		print( 'Out file:', os.path.basename(outFileStr))
	
	# check FASTA
	fastaIndex = checkFasta( fastaFileStr )
	
	# check BED
	valid = checkBed( bedFileStr )
	if not valid:
		print('ERROR: Incorrect BED file. BED must have score field.')
		exit()
	
	fastaFile = pysam.Fastafile( fastaFileStr )
	dmrFile = pybedtools.BedTool( bedFileStr )
	
	# build position dict if necessary
	if isPos:
		if isPrint:
			print('Getting positions')
		posDict = readPosFile( posFileStr )
	else:
		posDict = None
	
	#updatedDmrs = dmrFile.each( computeCGDens, fasta=fastaFile )
	if isPrint:
		print( dmrFile.count(), 'features in input file' )
	updatedDmrs = dmrFile
	
	if isPrint:
		print('Running process')
		
	nTries, lastFile, leftOverDMRs = runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, posDict, posThresh, bamFileAr, minReads, seedVal, includeFileStr, isPrint )
	if nTries >= maxTries:
		print('WARNING: Did not converge')
		leftOverDMRs.saveas( outId + '_unconverged.bed')
	
	# save final file
	
	pybedtools.BedTool(lastFile).sort(faidx=fastaIndex).moveto(outFileStr)
	
	if not isKeep:
		tmpFileAr = glob.glob('tmp*')
		for x in tmpFileAr:
			os.remove(x)
	if isPrint:
		print( 'Done.' )

def readBamList( bamListFile ):
	inFile = open( bamListFile, 'r' )
	outAr = []
	
	for line in inFile:
		if line.startswith('#'):
			continue
		bamFileStr = line.rstrip()
		bamIndex = bamFileStr + '.bai'
		if os.path.isfile(bamFileStr) == False:
			print( 'WARNING: BAM file {:s} does not exist...skipping'.format( os.path.basename(bamFileStr) ) )
			continue
		if os.path.isfile(bamIndex) == False:
			print( 'Creating BAM index for', os.path.basename(bamFileStr))
			pysam.index( bamFileStr )
		outAr += [bamFileStr]
	# end for line
	inFile.close()
	return outAr
	
def checkFasta( fastaFileStr ):
	fastaIndex = fastaFileStr + '.fai'
	if os.path.isfile(fastaFileStr) == False:
		print('ERROR: FASTA file does not exist')
		exit()
	elif os.path.isfile(fastaIndex) == False:
		print('WARNING: FASTA index file does not exist...creating')
		pysam.faidx( fastaFileStr )
	return fastaIndex
	
def checkBed( bedFileStr ):
	with open( bedFileStr ) as f:
		first_line = f.readline()
		lineAr = first_line.rstrip().split('\t')
		return len(lineAr) >= 5

def readPosFile( posFileStr ):
	'''
		Pos file is currently 1-indexed and we need to convert it to 0-indexed
	'''
	if os.path.isfile(posFileStr) == False:
		print('ERROR: position file does not exist')
		exit()
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
		#outDict[chrm] += [pos]
		bisect.insort( outDict[chrm], pos )
	# end for line
	posFile.close()
	return outDict

def runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, posDict, posThresh, bamFileAr, minReads, seedVal, includeFileStr, isPrint ):
	
	excludeFileStr = 'tmp.exclude.bed'
	outFileStr = 'tmp.equiv_-1.bed'
	remainFileStr = 'tmp.original.bed'
	remainingDMRs = updatedDmrs.saveas(remainFileStr)
	
	foundRegions = pybedtools.BedTool('', from_string=True).saveas(outFileStr)
	
	i=0
	while i < maxTries:
		x = foundRegions.cat('tmp.original.bed', force_truncate=True).moveto(excludeFileStr)
		
		# try to find intervals
		tryAgainDMRs, newRegions = runTry(i, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, posDict, posThresh, bamFileAr, minReads, includeFileStr, isPrint )
		
		#  add new found regions to old
		
		foundRegions = pybedtools.BedTool(newRegions).cat(outFileStr, postmerge=False)
		
		outFileStr = 'tmp.equiv_{:d}.bed'.format(i)
		remainFileStr = 'tmp.remain_{:d}.bed'.format(i)
		foundRegions.saveas(outFileStr)
		
		if len(tryAgainDMRs) == 0: # none left
			break
		
		remainingDMRs = remainingDMRs.filter(lambda x: x.name in tryAgainDMRs).saveas(remainFileStr)
		i+= 1
		
	# end for i
	return i, outFileStr, remainingDMRs
	
def runTry( i, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, posDict, posThresh, bamFileAr, minReads, includeFileStr, isPrint ):
	
	isPos = posThresh != -1
	isBam = len(bamFileAr) > 0
	
	if isBam:
		bamPointerAr = [ pysam.AlignmentFile( bamFileStr, 'rb' ) for bamFileStr in bamFileAr ]
		
	remainDMRs = pybedtools.BedTool(remainFileStr)

	# generate shuffle of remaining DMRs
	if seedVal == None and includeFileStr == None:
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr)
	elif seedVal == None:
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, incl=includeFileStr )
	elif includeFileStr == None:
		seedVal += i
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, seed=seedVal )
	else:
		seedVal += i
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, incl=includeFileStr, seed=seedVal )
	
	# iterate to see which to keep
	tryAgainDMR = [] # DMR labels
	outAr = [] # output regions
	
	for feature in outBed:
		nRead = 0
		isValid = True
		# if close enough
		#print('  ', cgDens, feature.score, abs(cgDens - float(feature.score) ) <= THRESH)
		if not isPos and not isBam: # no thresholds
			outAr += [ feature ]
			continue

		if isPos: # also check position
			posAr = posDict.get( feature.chrom )
			nPos = _computePos( posAr, feature.start, feature.stop )
			#print(feature.name, nPos)
			isValid = isValid and (nPos >= posThresh)
		
		if isBam: # check read coverage
			nReadAr = [ bamFile.count(feature.chrom, feature.start, feature.stop) for bamFile in bamPointerAr ]
			nRead = min ( nReadAr )
			isValid = isValid and (nRead >= minReads)
		
		if isValid:
			#acceptDMR += [feature.label]
			feature.score = nRead
			outAr += [ feature ]
		else:
			 tryAgainDMR += [ feature.name ]
	# end for feature
	# close bam files
	if isBam:
		t = [ x.close() for x in bamPointerAr ]
		
	if isPrint:
		print('{:d} regions remaining'.format( len(tryAgainDMR)) )
	return tryAgainDMR, outAr

def _computePos( posAr, start, end ):
	if posAr == None: # no positions on chrm
		return 0
	# index of first position >= start
	sIndex = bisect_ge( posAr, start )
	# index of last position < end
	eIndex = bisect_lt( posAr, end )
	
	if eIndex == None or sIndex == None or eIndex <= sIndex:
		return 0
	# how many positions are there
	return eIndex - sIndex

'''def readCoverage( bamPointerAr, chrm, start, end ):
	outAr = []
	
	for bamFile in bamPointerAr:
		nReads = bamFile.count( chrm, start, end )'''
		
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
	
def parseInputs( argv ):
	maxTries = NTRIES
	posThresh = POSTHRESH
	minReads = NREADS
	bamFileList = None
	outId = None
	seedVal = None
	posFileStr = None
	includeFileStr = None
	isPrint = True
	isKeep = False
	startInd = 0
	
	for i in range(min(9, len(argv))):
		if argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-k':
			isKeep = True
			startInd += 1
		elif argv[i].startswith('-n='):
			try:
				maxTries = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: max tries must be an integer...using default', NTRIES)
			startInd += 1
		elif argv[i].startswith('-s='):
			try:
				seedVal = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: seed value must be an integer...using None')
			startInd += 1
		elif argv[i].startswith('-t='):
			posFileStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-p='):
			try:
				posThresh = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: position threshold must be an integer...using None' )
			startInd += 1
		elif argv[i].startswith('-r='):
			try:
				minReads = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: min read coverage must be an interger...using default', NREADS)
			startInd += 1
		elif argv[i].startswith('-o='):
			outId = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-b='):
			bamFileList = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-i='):
			includeFileStr = argv[i][3:]
			startInd += 1
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith('-'):
			print('ERROR: invalid parameter', argv[i], 'use -h to see valid parameters')
			exit()
	# end for
	
	if (posFileStr == None and posThresh != -1) or (posFileStr != None and posThresh == -1):
		print('WARNING: need to specificy both pos file and pos threshold for position filtering to work...no position filtering')
		posThresh == -1
		posFileStr = None
	
	fastaFileStr = argv[startInd]
	bedFileStr = argv[startInd + 1]
	
	processInputs( bedFileStr, fastaFileStr, posFileStr, bamFileList, minReads, maxTries, outId, seedVal, posThresh, includeFileStr, isKeep, isPrint )

def printHelp():
	print( 'python\tgenome_bed_shuffle_v3.py [-h] [-q] [-k] [-s=seed_val]' )
	print( '\t[-n=max_tries] [-o=out_id] [-t=pos_file] [-p=pos_threshold] [-b=bam_files]')
	print('\t[-r=min_read_cov] [-i=include_file]<fasta_file> <bed_file>')
	print('WARNING: temporary files are written using predefined names. DO NOT run multiple instances of this script in the same directory at the same time.')

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
