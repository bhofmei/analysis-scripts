import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam

# Usage: region_read_allele.py [-q] [-h] [-z] [-o=out_id] <pos_file> <bed_file>  <bam_file>

#TRANSITIONS = ['uu', 'uM', 'Mu', 'MM','u.','M.','.u','.M']
TRANSITIONS = ['uu', 'uM', 'Mu', 'MM']

def processInputs( bamFileStr, bedFileStr, posFileStr, outID, isCompress, isPrint ):
	
	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'BAM file:', os.path.basename(bamFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
	
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
	for dmr in dmrPosAr:
		dmrRes = processRegion( bamFile, dmr, isPrint )
		outMatrix += [dmrRes]
	
	# close bam and fasta objects
	bamFile.close()
	
	# write output
	if outID == None:
		outID = os.path.basename(bamFileStr).replace('.bam', '')
	outFileStr = outID + '_dmr_reads.tsv'
	outFileStr += '.gz' if isCompress else ''
	if isPrint:
		print( 'Writing output to', outFileStr)
	
	info = '#from_script: region_haplo_allele.py; bed_file:{:s}; starting_regions:{:d}; remaining_regions:{:d}\n'.format(os.path.basename(bedFileStr), len(dmrAr), len(dmrPosAr))
	
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
			
def processRegion( bamFile, dmrRegion, isPrint ):
	chrm, start, end, label, posPlus, posMinus = dmrRegion
	outDict = {'plus': [], 'minus': []}
	
	# get reads
	reads = bamFile.fetch(chrm, start, end)
	
	# loop through reads
	for read in reads:
		rPos = read.get_reference_positions(full_length=True)
		rSeq = read.query_sequence
		isRev = read.is_reverse
		
		if isRev:
			a = encodeReverseRead( posMinus, rPos, rSeq )
			outDict['minus'] += ( [a] if a != None else [] )
		else:
			a = encodeForwardRead( posPlus, rPos, rSeq )
			outDict['plus'] += ( [a] if a != None else [] )
	map = buildMap(len(posPlus), outDict['plus'])
	map2 = buildMap(len(posMinus), outDict['minus'])
	#n = max(len(posPlus), len(posMinus))
	#minusFirst = False
	#if(posMinus[0]<posPlus[0]):
	#	minusFirst = True
		#print('warning adjust map -', label, posPlus[0], posMinus[0])
	#if (abs(len(posPlus) - len(posMinus)) > 1):
	#	print('warning', label, len(posPlus), len(posMinus))
	#map3 = buildMap2( n, outDict['plus'], outDict['minus'], minusFirst)
	
	# find paths
	ppaths = findPaths( map )
	mpaths = findPaths( map2 )
	#fPaths = list( filter(lambda t: len(t[1]) == n, paths ) )
	pospaths = list( filter(lambda t: len(t[1]) == len(posPlus), ppaths ) )
	minpaths = list( filter(lambda t: len(t[1]) == len(posMinus), mpaths ) )
	
	pospaths.sort( reverse = True )
	minpaths.sort( reverse = True )
	# end for read
	outStr = '\n### {:s}\t\t{:s} : {:d} - {:d}\n'.format(label, chrm, start, end)
	outStr += '+ ' + str(posPlus) + '\n- ' + str(posMinus) + '\n'
	for x in outDict['plus']:
		outStr += '+|' + x + '|\n'
	for y in outDict['minus']:
		outStr += '-|' + y + '|\n'
	outStr += '# + #\n'
	for x in map:
		outStr += str(x) + '\n'
	for i in range(min(len(pospaths), 10)):
		tup = pospaths[i]
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:s}\n'.format(chrm, start, end, label, i, tup[0], tup[1])
	#outStr += '\n'
	outStr += '# - #\n'
	for y in map2:
		outStr += str(y) + '\n'
	for i in range(min(len(minpaths), 10)):
		tup = minpaths[i]
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:s}\n'.format(chrm, start, end, label, i, tup[0], tup[1])
	outStr += '\n'
	#for z in map3:
	#	outStr += str(z) + '\n'
	#outStr += '\n'
	#outStr = ''
	#for i in range(min(len(fPaths), 10)):
	#	tup = fPaths[i]
	#	outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:s}\n'.format(chrm, start, end, label, i, tup[0], tup[1])
	
	return outStr	

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

def buildMap(n, plusReads):
	
	counts = [[0]*len(TRANSITIONS) for x in range(n-1)]
	
	for read in plusReads:
		for j in range(n-1):
			ch = read[j:j+2]
			#print(ch)
			try:
				i = TRANSITIONS.index(ch)
				#print(i)
				counts[j][i]+= 1
			except ValueError:
				pass
		# end for j
	# end for read
	return counts

def buildMap2( n, plusReads, minusReads, minusFirst ):
	
	counts = [[0]*len(TRANSITIONS) for x in range(n-1)]
	plusStart = (1 if minusFirst else 0)
	minusStart = (0 if minusFirst else 1)
	for read in plusReads:
		for j in range(plusStart, n-1):
			ch = read[j:j+2]
			#print(ch)
			try:
				i = TRANSITIONS.index(ch)
				#print(i)
				counts[j][i]+= 1
			except ValueError:
				pass
		# end for j
	# end for read
	for read in minusReads:
		for j in range(minusStart, n-1):
			ch = read[j:j+2]
			try:
				i = TRANSITIONS.index(ch)
				counts[j][i] += 1
			except ValueError:
				pass
		# end for j
	# end for read
	return counts

def findPaths(matrix):
	outAr = []
	if len(matrix) == 0:
		return outAr
	for i in range(len(TRANSITIONS)):
		x = matrix[0][i]
		if x != 0:
			outAr += pathHelper(matrix[1:], TRANSITIONS[i], x)
	return outAr

def pathHelper( matrix, path, count ):
	# base case -> end of matrix
	if len(matrix) == 0:
		return [(count, path)]
	# otherwise go down the branches
	n = path[-1]
	curRow = matrix[0]
	if n == 'u':
		#testPaths = [0, 1, 4]
		testPaths = [0, 1]
	elif n == 'M':
		#testPaths = [2, 3, 5]
		testPaths = [2, 3]
	#elif n == '.':
		#testPaths = [6, 7]
	else:
		testPaths = []
	outAr = []
	for p in testPaths:
		if curRow[p] != 0:
			nextPath = path[:-1] + TRANSITIONS[p]
			newCount = count + curRow[p]
			newMatrix = ([] if len(matrix)== 1 else matrix[1:])
			m = pathHelper(matrix[1:], nextPath, newCount)
			outAr += m
	return outAr
	
def writeOutput( outFileStr, outMat, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	outFile.write(info)
	headerAr = ['chrm', 'start', 'end', 'label', 'i', 'count', 'seq']
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
	bamFile = argv[startInd +2]
	
	
	processInputs( bamFile, bedFile, posFile, outID, isCompress, isPrint )
	
def printHelp():
	print('Usage: python3 region_haplo_allele.py [-q] [-h] [-o=out_id] <pos_file> <bed_file> <bam_file>')


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
