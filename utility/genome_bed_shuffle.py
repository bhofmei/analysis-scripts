import sys, math, glob, multiprocessing, subprocess, os, bisect, random, re
import pysam, pybedtools

# Usage: python genome_bed_shuffle.py [-h] [-q] [-t=gc_threshold] [-s=seed_val] [-n=max_tries] [-o=out_id] <fasta_file> <bed_file>
NTRIES=1000
THRESH=0.05

def processInputs( fastaFileStr, bedFileStr, maxTries, gcThresh, outId, seedVal, isPrint ):
	
	if outId == None:
		tmp = os.path.basename( bedFileStr )
		outId = re.sub( "\.bed", "", tmp, flags=re.IGNORECASE )
	outFileStr = outId + '_equiv' + '.bed'
		
	if isPrint:
		print( 'Max Tries:', maxTries )
		print( 'GC Threshold:', gcThresh )
		print( 'Starting seed:', seedVal )
		print( 'Fasta file:', os.path.basename(fastaFileStr) )
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
	
	if isPrint:
		print( 'Getting regions and nucleotide density' )
	
	updatedDmrs = dmrFile.each( computeCGDens, fasta=fastaFile )
	
	nTries, lastFile, leftOverDMRs = runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, gcThresh, seedVal, isPrint )
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

def runProcess( maxTries, updatedDmrs, fastaFile, fastaIndex, gcThresh, seedVal, isPrint ):
	
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
		tryAgainDMRs, newRegions = runTry(i, gcThresh, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, isPrint )
		
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
	
def runTry( i, gcThresh, seedVal, remainFileStr, excludeFileStr, fastaFile, fastaIndex, isPrint ):
	
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
		# new CG Density
		cgDens = _computeCG( fastaFile, feature.chrom, feature.start, feature.stop )
		# if close enough
		#print('  ', cgDens, feature.score, abs(cgDens - float(feature.score) ) <= THRESH)
		if (gcThresh==-1) or (abs(cgDens - float(feature.score) ) <= gcThresh):
			#acceptDMR += [feature.label]
			outAr += [ feature ]
		else:
			 tryAgainDMR += [ feature.name ]
	# end for feature

	if isPrint:
		print('{:d} regions remaining'.format( len(tryAgainDMR)) )
	return tryAgainDMR, outAr

def parseInputs( argv ):
	maxTries = NTRIES
	gcThresh = THRESH
	outId = None
	seedVal = None
	isPrint = True
	startInd = 0
	
	for i in range(min(5, len(argv))):
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
	fastaFileStr = argv[startInd]
	bedFileStr = argv[startInd + 1]
	
	processInputs( fastaFileStr, bedFileStr, maxTries, gcThresh, outId, seedVal, isPrint )

def printHelp():
	print( 'Usage: python genome_bed_shuffle.py [-h] [-q] [-t=gc_threshold] [-s=seed_val] [-n=max_tries] [-o=out_id] <fasta_file> <bed_file>' )
	print('WARNING: temporary files are written using predefined names. DO NOT run multiple instances of this script in the same directory at the same time.')

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
