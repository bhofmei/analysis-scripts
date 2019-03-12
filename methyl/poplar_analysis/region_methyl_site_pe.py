import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
# import pysam, pybedtools

# Usage: python3 region_methyl_site_pe.py [-h] [-q] [-z] [-m] [-o=out_id] <pos_file> <bed_file> <allc_file>

NUMPROC = 1

def processInputs( allcFileStr, bedFileStr, posFileStr, outID, numProc, isMerge, isCompress, isPrint ):
	if isPrint:
		print( 'Position file:', os.path.basename(posFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
		print( 'Allc file:', os.path.basename(allcFileStr) )
	
	# get DMRs
	dmrAr = readBED( bedFileStr )
	if isPrint:
		print( len(dmrAr), 'regions found in', os.path.basename(bedFileStr) )
	
	# get keep positions - 0-based
	if isPrint:
		print( 'Reading position file', os.path.basename(posFileStr) )
	posDict = readPosFile( posFileStr )
	
	# get allC information
	if isPrint:
		print( 'Reading allc file', os.path.basename(allcFileStr) )
	allcDict = readAllcPos( allcFileStr, posDict )
	
	# determine keep positions in each DMR
	if isPrint:
		print( 'Intersecting regions and positions' )
	dmrPosAr = intersectRegionsPos( dmrAr, posDict, isMerge )
	if isPrint:
		print( 'Processing remaining', len(dmrPosAr), 'regions' )
	
	outMatrix = []
	for dmr in dmrPosAr:
		chrm = dmr[0]
		tmpDict = allcDict.get(chrm)
		if tmpDict != None:
			outStr = processRegion( dmr, tmpDict, isMerge )
			outMatrix += [ outStr ]
		# end if
	# end for
	
	if outID == None:
		outID = os.path.basename(allcFileStr).replace('.tsv', '').replace('.gz','')
		
	outFileStr = outID + '_dmr_methyl-site.tsv'
	outFileStr += '.gz' if isCompress else ''
	
	info = '#from_script: region_methyl_site.py; allc_file:{:s}; bed_file:{:s}; remaining_regions:{:d}; merge_strands:{:s}\n'.format(os.path.basename(allcFileStr), os.path.basename(bedFileStr), len(dmrPosAr), str(isMerge))
	
	if isPrint:
		print( 'Writing output to', outFileStr)
	writeOutput( outFileStr, outMatrix, info, isCompress )
	
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

def readAllcPos( allcFileStr, posDict ):
	outDict = {}
	
	allcFile = openAllc( allcFileStr )
	currChrm = None
	currPosAr = None
	
	for line in allcFile:
		lineAr = line.rstrip().split('\t')
		# header
		if line.startswith('#') or len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		# (0) chrm (1) pos 1-index (2) strand (3) mc class (4) mc count (5) tc count
		# (6) isMethyl
		chrm = lineAr[0]
		pos = int( lineAr[1] ) - 1
		mc = int( lineAr[4] )
		tc = int( lineAr[5] )
		# we've moved to next chrm
		if currChrm != chrm:
			currPosAr = posDict.get(chrm)
			outDict[chrm] = {}
			currChrm = chrm
		
		# see if this pos is in currPosAr
		i = (None if currPosAr == None else bisect_index(currPosAr, (pos, '')))
		if i != None:
			# get data and add to output
			outDict[chrm][pos] = (mc, tc)
			
		# i == None, not in list, ignore it and move on
	# end for line
	allcFile.close()
	return outDict

def openAllc( allcFileStr ):
	if allcFileStr.endswith('.gz'):
		return gzip.open( allcFileStr, 'rt' )
	else:
		return open( allcFileStr, 'r' )

def processChrm( dmrAr, allcDict, isMerge ):
	outStr = ''
	
	for dmr in dmrAr:
		tmpStr = processRegion( dmr, allcDict, isMerge )
		outStr += tmpStr
	# end for
	
	return outStr

def intersectRegionsPos( dmrAr, posDict, isMerge ):
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
		
		# check start pos is + and end is -, otherwise adjust sIndex and eIndex
		if isMerge:
			o = posAr[sIndex]
			if posAr[sIndex][1] == '-' and posAr[sIndex-1][1] == '+':
				sIndex -= 1 # start with the one before
			elif posAr[sIndex][1] == '-':
				print('Warning: error with symmetry', posAr[sIndex])
			p = posAr[eIndex]
			if posAr[eIndex][1] == '+' and posAr[eIndex+1][1] == '-':
				eIndex += 1 # end with one after
			elif posAr[eIndex][1] == '+':
				print('Warning: error with symmetry', posAr[eIndex])
			#print(o, '->', posAr[sIndex], '--', p, '->', posAr[eIndex])
		# get the positions
		cPos = posAr[ sIndex:eIndex+1 ]
		#print('c', cPos)
		
		# not enough positions, move on
		if len(cPos) <= 2:
			continue
		
		# not merge, return just the positions
		if isMerge:
			n = len(cPos)
			if n % 2 != 0:
				print('Warning: uneven position list', cPos)
				outPos = []
			else:
				n1 = int(n / 2)
				outPos = [(cPos[i*2][0], cPos[i*2+1][0]) for i in range(n1)]
		else:
			outPos = [t[0] for t in cPos]
		
		if outPos != []:
			outAr += [(chrm, start, end, label, outPos)]
	
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
		if i != len(a) and a[i][0] == x[0]:
			return i
		return None
	except TypeError:
		if x in a:
			return a.index(x)
		else:
			return None



def processRegion( dmr, allcDict, isMerge ):
	chrm, start, end, label, posAr = dmr
	outStr = ''
	for i, pos in enumerate(posAr):
		if isMerge:
			mc1, tc1 = allcDict[pos[0]]
			mc2, tc2 = allcDict[pos[1]]
			mc = mc1 + mc2
			tc = tc1 + tc2
			p = pos[0]
		else:
			# not merged
			mc, tc = allcDict[pos]
			p = pos
		
		wm = float(mc) / float(tc)
		outStr += '{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:f}\n'.format( chrm, start, end, label, i, p, mc, tc, wm)
	# end for pos
	
	return outStr	

def writeOutput( outFileStr, outMatrix, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	headerAr = ['chrm', 'start', 'end', 'label', 'i', 'pos', 'mc', 'tc', 'wm']
	outFile.write( info )
	outFile.write( '\t'.join(headerAr) + '\n' )
	
	for x in outMatrix:
		outFile.write( x )
	# end for x
	outFile.close()

def parseInputs( argv ):
	outID = None
	numProc = NUMPROC
	isMerge = False
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
		elif argv[i] == '-m':
			isMerge = True
			startInd += 1
		elif argv[i].startswith('-o='):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'WARNING: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
			startInd += 1
		elif argv[i].startswith('-'):
			print('ERROR: "{:s}" is not a valid parameter\nUse -h to see valid parameters'.format(argv[i]))
	# end for
	
	posFileStr = argv[startInd]
	bedFileStr = argv[startInd+1]
	allcFileStr = argv[startInd +2]
	
	processInputs( allcFileStr, bedFileStr, posFileStr, outID, numProc, isMerge, isCompress, isPrint )

def printHelp():
	print('Usage: python3 region_methyl_site.py [-h] [-q] [-z] [-m] [-o=out_id] [-p=num_proc] <pos_file> <bed_file> <allc_file>')
	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
