import sys, math, glob, multiprocessing, subprocess, os, bisect, random, re
import pysam, pybedtools

# Usage: python genome_bed_shuffle.py [-h] [-q] [-t=gc_threshold] [-s=seed_val] [-n=max_tries] [-o=out_id] [-p=pos_file] [-r=pos_threshold] [-i=include_file] <fasta_file> <bed_file>
NTRIES=1000
THRESH=0.05
POSTHRESH = -1

def processInputs( bedFileStr, fastaFileStr, posFileStr, includeFileStr, maxTries, gcThresh, outId, seedVal, posThresh, isPrint ):

	isPos = (posFileStr != None)
	isInc = (includeFileStr != None)
	
	if outId == None:
		tmp = os.path.basename( bedFileStr )
		outId = re.sub( "\.bed", "", tmp, flags=re.IGNORECASE )
	outFileStr = outId + '_equiv' + '.bed'
		
	if isPrint:
		print( 'Max Tries:', maxTries )
		print( 'GC Threshold:', gcThresh )
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
		if isPos:
			print( 'Position Threshold:', posThresh)
			print( 'Position file:', os.path.basename(posFileStr) )
		if isInc:
			print( 'Include file:', os.path.basename(includeFileStr) )
		print( 'Starting seed:', seedVal )
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
		posDict = readPosFile( posFileStr )
	else:
		posDict = None
	
	# check include bed file exist if necessary
	if isInc and os.path.isfile(includeFileStr) == False:
		print('ERROR: Include BED file is specificed but does not exist')
		exit()
	
	if isPrint:
		print( 'Getting regions and nucleotide density' )
	
	updatedDmrs = dmrFile.each( computeCGDens, fasta=fastaFile )
	
	if isPrint:
		print('Running process')
		
	nTries, lastFile, leftOverDMRs = runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, includegcThresh, posDict, posThresh, seedVal, isPrint )
	if nTries >= maxTries:
		print('WARNING: Did not converge')
		leftOverDMRs.saveas( outId + '_unconverged.bed')
	
	# save final file
	
	pybedtools.BedTool(lastFile).sort(faidx=fastaIndex).moveto(outFileStr)
	
	tmpFileAr = glob.glob('tmp*')
	#print(tmpFileAr)
	for x in tmpFileAr:
		os.remove(x)
	if isPrint:
		print( 'Done.' )
	
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
	
def computeCGDens( feature, fasta=None ):
	if fasta == None:
		score = -1
	else:
		score = _computeCG( fasta, feature.chrom, feature.start, feature.stop )
	feature.score = str(score)
	return feature

def _computeCG( fastaFile, chrm, start, end ):
	refSeq = fastaFile.fetch( chrm, start, end)
	
	cgList = list(filter(lambda r: r == 'C' or r == 'G', refSeq))
	return len(cgList) / float(end - start)

def runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, gcThresh, posDict, posThresh, seedVal, isPrint ):
	
	excludeFileStr = 'tmp.exclude.bed'
	outFileStr = 'tmp.equiv_-1.bed'
	remainFileStr = 'tmp.original.bed'
	remainingDMRs = updatedDmrs.saveas(remainFileStr)
	#outBedFile = 'equiv_regions.bed'
	
	foundRegions = pybedtools.BedTool('', from_string=True).saveas(outFileStr)
	
	i=0
	while i < maxTries:
		x = foundRegions.cat('tmp.original.bed', force_truncate=True).moveto(excludeFileStr)
		
		# try to find intervals
		tryAgainDMRs, newRegions = runTry(i, gcThresh, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, posDict, posThresh, isPrint )
		
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
	
def runTry( i, gcThresh, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, posDict, posThresh, isPrint ):
	
	isGC = gcThresh != -1
	isPos = posThresh != -1
	remainDMRs = pybedtools.BedTool(remainFileStr)
	#if isPrint:
		#print( 'Iteration {:d} with {:d} regions...'.format(i, remainDMRs.count()), end='' )
	# generate shuffle of remaining DMRs
	if seedVal == None:
		#outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, chrom=True )
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr)
	else:
		#outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, chrom=True, seed=seedVal )
		seedVal += i
		outBed = remainDMRs.shuffle( g=fastaIndex, excl=excludeFileStr, seed=seedVal )
	
	# iterate to see which to keep
	#acceptDMR = [] # DMR labels
	tryAgainDMR = [] # DMR labels
	outAr = [] # output regions
	
	for feature in outBed:
		isValid = True
		# if close enough
		#print('  ', cgDens, feature.score, abs(cgDens - float(feature.score) ) <= THRESH)
		if not isGC and not isPos: # no thresholds
			outAr += [ feature ]
			continue
		
		if isGC: # check gc content if necessary
			cgDens = _computeCG( fastaFile, feature.chrom, feature.start, feature.stop )
			isValid = (abs(cgDens - float(feature.score) ) <= gcThresh)
		
		if isPos: # also check position
			posAr = posDict.get( feature.chrom )
			nPos = _computePos( posAr, feature.start, feature.stop )
			#print(feature.name, nPos)
			isValid = isValid and (nPos >= posThresh)
		
		if isValid:
			#acceptDMR += [feature.label]
			outAr += [ feature ]
		else:
			 tryAgainDMR += [ feature.name ]
	# end for feature

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
	gcThresh = THRESH
	posThresh = POSTHRESH
	outId = None
	seedVal = None
	posFileStr = None
	isPrint = True
	startInd = 0
	
	for i in range(min(7, len(argv))):
		if argv[i] == '-q':
			isPrint = False
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
		elif argv[i].startswith('-p='):
			posFileStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-r='):
			try:
				posThresh = int(argv[i][3:])
			except ValueError:
				print( 'WARNING: position threshold must be an integer...using None' )
			startInd += 1
		elif argv[i].startswith('-t='):
			try:
				gcThresh = float(argv[i][3:])
				if gcThresh >= 1:
					gcThresh = gcThresh / 100.0
			except ValueError:
				print( 'WARNING: GC threshould must be a float...using default', THRESH)
			startInd += 1
		elif argv[i].startswith('-o='):
			outId = argv[i][3:]
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
	
	processInputs( bedFileStr, fastaFileStr, posFileStr, maxTries, gcThresh, outId, seedVal, posThresh, isPrint )

def printHelp():
	print( 'Usage: python genome_bed_shuffle.py [-h] [-q] [-t=gc_threshold] [-s=seed_val] [-n=max_tries] [-o=out_id] [-p=pos_file] [-r=pos_threshold] <fasta_file> <bed_file>' )
	print('WARNING: temporary files are written using predefined names. DO NOT run multiple instances of this script in the same directory at the same time.')

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
