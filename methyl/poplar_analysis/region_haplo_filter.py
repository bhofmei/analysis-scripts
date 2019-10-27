import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam
from DmrPos import *
from bisect_func import *

# Usage: python3 region_haplo_filter.py [-q] [-z] [-h] [-r=min_reads] [-p=min_pos]
# [-o=out_id] <pos_file> <bed_file> <bam_file> [bam_files*]

MINREADS=10
MINPOS=2

def processInputs( bamFileAr, posFileStr, bedFileStr, outId, minReads, minPos, isCompress, isPrint ):

	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'BAM files:', ','.join([os.path.basename(x) for x in bamFileAr] ))
		print( 'BED file:', os.path.basename(bedFileStr) )
		print( 'Min reads: {:d}'.format(minReads) )
		print( 'Min position overlap per strand: {:d}'.format(minPos) )
	
	# 1. get DMRs
	dmrAr = readBED( bedFileStr )
	if isPrint:
		print( len(dmrAr), 'regions found')
		
	# 2. get positions
	posDict, nPos = readPosFile( posFileStr )
	if isPrint:
		print( nPos, 'positions found' )
		
	# 3. filter by num positions
	if isPrint:
		print( 'Filtering regions by min number of positions...' )
	dmrPosAr = filterRegionsPos( dmrAr, posDict, minPos ) 
	if isPrint:
		print( len(dmrPosAr), 'regions remaining' )
	
	# dmrPosAr tuples have: (chrm, start, end, label, posPlusAr, posMinusAr)
	
	# 4. check bams
	for bamFileStr in bamFileAr:
		checkBAM( bamFileStr )
	
	# 5. filter by min reads
	if isPrint:
		print( 'Filtering regions by min number of reads in all samples...' )
	for bamFileStr in bamFileAr:
		dmrPosAr = filterRegionsBam( dmrPosAr, bamFileStr, minReads )
	if isPrint:
		print( len(dmrPosAr), 'regions remaining' )
	
	# 6. report results
	if outId == None:
		outId = os.path.basename(bedFileStr).replace('.bed', '').replace('.BED','')
	
	outFileStr = outId + '_haplo_regions.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr )
	
	info = '#from_script: region_haplo_filter.py; bed_file:{:s}; pos_file:{:s}; min_pos:{:d}; min_reads:{:d}; starting_regions:{:d}; remaining_regions:{:d}'.format( os.path.basename(bedFileStr), os.path.basename(posFileStr), minPos, minReads, len(dmrAr), len(dmrPosAr))
	
	writeOutput( outFileStr, dmrPosAr, info, isCompress)
	if isPrint:
		print('Done')

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
		#outAr += [(chrm, start, end, dmrName)]
		t = DmrPos(chrm, start, end, dmrName)
		outAr += [t]
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
	nPos = 0
	for line in posFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or lineAr[1].isdigit() == False:
			continue
		chrm = lineAr[0]
		pos = int(lineAr[1]) - 1
		strand = lineAr[2]
		
		if outDict.get(chrm) == None:
			outDict[chrm] = []
		bisect.insort( outDict[chrm], (pos,strand) )
		nPos += 1
		
	# end for line
	posFile.close()
	return outDict, nPos
	
def filterRegionsPos( dmrAr, posDict, minPos ):
	outAr = []
	
	# loop through regions
	for region in dmrAr:
		#chrm, start, end, label = region
		chrm, start, end = region.getCoord()
		
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
		# if too few positions per strand, don't keep it
		if len(posPlus) < minPos or len(posMinus) < minPos:
			continue
		#outAr += [(chrm, start, end, label, posPlus, posMinus)]
		region.setPos( posPlus, posMinus )
		outAr += [region]
	# end for region
	return outAr
	
def checkBAM( bamFileStr ):
	bamIndex = bamFileStr + '.bai'
	if os.path.isfile(bamFileStr) == False:
		print('ERROR: BAM file does not exist')
		exit()
	elif os.path.isfile(bamIndex) == False:
		print('WARNING: BAM index file does not exist...creating')
		pysam.index( bamFileStr )
	return True

def filterRegionsBam( dmrPosAr, bamFileStr, minReads ):
	outAr = []
	
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	for region in dmrPosAr:
		chrm, start, end = region.getCoord()
		nReads = bamFile.count( chrm, start, end )
		if nReads >= minReads:
			outAr += [region]
	# end for region
	bamFile.close()
	return outAr

def writeOutput( outFileStr, dmrPosAr, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info + '\n')
	headerAr = ['chrm', 'start', 'end', 'label', 'posAr', 'negAr']
	outFile.write('\t'.join(headerAr) + '\n' )
	
	dmrPosAr.sort()
	
	for region in dmrPosAr:
		outFile.write(str(region) + '\n')
	# end for region
	outFile.close()
	
def parseInputs( argv ):
	outId = None
	isPrint = True
	isCompress = False
	minReads = MINREADS
	minPos = MINPOS
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
			outId = argv[i][3:]
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
	# end for i
	
	posFileStr = argv[startInd]
	bedFileStr = argv[startInd+1]
	
	bamFileAr = []
	for j in range(startInd+2, len(argv)):
		bamFileAr += [argv[j]]
	
	processInputs( bamFileAr, posFileStr, bedFileStr, outId, minReads, minPos, isCompress, isPrint )

def printHelp():
	print('Usage:\tpython3 region_haplo_filter.py [-q] [-h] [-r=min_reads] [-p=min_pos]')
	print( '\t[-o=out_id] <pos_file> <bed_file> <bam_file> [bam_files*]' )
	print()
	print('Note: include ALL bam files to be used')
	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
