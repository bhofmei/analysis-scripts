import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 determine_cg_smps.py [-o=outPre][-c=chromosome_names] <min_cov> <allC_path> <line1_name> <line2_name> [lineN_name]*
# Using allC files, determines the CG-SMPs for all samples (outputs these positions) and computes a distance matrix between all samples based on these positions
# CG smp defined as: 1) methylation not the same among all samples and 2) all samples have at least a) minCov number reads or b) methylated read at the position

CHRMNAMES=['Chr1','Chr2','Chr3','Chr4','Chr5']
#CHRMNAMES=['Chr1']
NUMPROC = 2

def processChrm( lineNamesAr, allCPath, chrm, minCov ):
	print( 'Analyzing chromosome {:s}...'.format(chrm) )
	chrmDict = generateChrmDict( lineNamesAr, allCPath, chrm )
	print( 'Cleaning dictionary for chromosome {:s}...'.format( chrm ) )
	chrmDict = cleanChrmDict( chrmDict, len(lineNamesAr), minCov)
	#print( chrmDict )
	print( 'Computing SMP difference matrix for chromosome {:s}...'.format( chrm ) )
	chrmMatrix = computeDifferencesChrm( chrmDict, len(lineNamesAr) )
	posAr = ['{:s}\t{:d}'.format( chrm, p ) for p in sorted(chrmDict.keys() ) ]
	return ( chrmMatrix, posAr )
	
def generateChrmDict( lineNamesAr, allCPath, chrm ):
	chrmDict = {}
	
	for line in lineNamesAr:
		allCFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, line, chrm ))
		chrmDict = readAllCtoDict( allCFileStr, chrmDict )
	return chrmDict

def cleanChrmDict( chrmDict, numLines, minCov ):
	
	# loop through each position (key)
	for pos in sorted( chrmDict.keys() ):
		dictEntry = chrmDict[pos]
		methAr = dictEntry[0]
		covAr = dictEntry[1]
		
		# check all lines had that position
		if len(methAr) != numLines:
			del chrmDict[pos]
			
		# check all the same -> remove
		# by seeing if all methylation states are equal to the first
		elif methAr.count( methAr[0] ) == len( methAr ):
			del chrmDict[pos]
		
		# check coverage
		# remove positions where a line had less reads than coverage and
		# position not methylated there
		else:
			lowCovMeth = [ methAr[i] for i in range(len(covAr)) if covAr[i] < minCov ]
			if lowCovMeth.count(0) > 0:
				del chrmDict[pos]
	return chrmDict	

def computeDifferencesChrm( chrmDict, numLines ):
	# difference matrix (note: only top diagonal will be filled)
	diffMatrix = [ [0]*numLines for x in range(numLines) ]
	
	# iterate through positions
	for pos in sorted( chrmDict.keys() ):
		methAr = chrmDict[pos][0]
		# iterate through lines pairwise
		for i in range( numLines ):
			for j in range( i+1, numLines ):
				diffMatrix[j][i] += ( methAr[i] != methAr[j] )
	return diffMatrix

def addToDifferenceMatrix( oldMatrix, currMatrix ):
	
	if oldMatrix == []:
		return currMatrix
	
	for i in range(len(oldMatrix)):
		for j in range(i+1,len(oldMatrix[0])):
			oldMatrix[i][j] += currMatrix[i][j]
	return oldMatrix
	
def readAllCtoDict( allCFileStr, chrmDict ):
	'''
		reads the AllC file and adds methylation information about chromosome
		positions to the dictionary
		if dictionary not empty and position wasn't previously in dictionary, 
		position is not added
		dictionary set up where key is position and value is 2D array where
		first array is lines methylation status and second array is lines
		coverage
	'''
	if chrmDict == {}:
		emptyDict = True
	else:
		emptyDict = False
	
	allCFile = open( allCFileStr, 'r' )
	print( 'Reading {:s}...'.format( allCFileStr ) )
	
	for line in allCFile:
		# header
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		
		# ignore non-CG positions
		if lineAr[3].startswith('CG') == False:
			continue
			
		pos = int(lineAr[1])
		cov = int( lineAr[5] )
		methylated = int( lineAr[6] )
		
		# fill dictionary for first file
		if emptyDict:
			chrmDict[pos] = [ [methylated], [cov] ]
		# position wasn't in dictionary -> continue
		elif chrmDict.get( pos ) == None:
			continue
		# in dictionary -> add information
		else:
			chrmDict[pos][0] += [methylated]
			chrmDict[pos][1] += [cov]
	allCFile.close()
	return chrmDict

def writeSMPOutput( outFileStr, posAr, minCov, lineNamesAr ):
	#print( posAr )
	outFile = open( outFileStr, 'w' )
	outFile.write( '#Lines included: {:s}\n'.format( ', '.join(lineNamesAr) ) )
	outFile.write( '#Minimum coverage: {:d}\n'.format( minCov ) )
	for pos in posAr:
		outFile.write( pos + '\n' )
	outFile.close()

def writeDiffMatrixOutput( outFileStr, lineNamesAr, diffMatrix ):

	numLines = len(lineNamesAr)
	for i in range(numLines):
		for j in range(i+1):
			diffMatrix[j][i] = -1
	
	outFile = open( outFileStr, 'w' )
	# header
	header = '\t' + '\t'.join( lineNamesAr )
	outFile.write( header + '\n' )
	
	for i in range(len(lineNamesAr) ):
		strAr = [ ('NA' if v == -1 else '{:d}'.format( v ) ) for v in diffMatrix[i] ]
		outStr = lineNamesAr[i].replace('-','.') + '\t' + '\t'.join( strAr ) + '\n'
		outFile.write( outStr )
	outFile.close()

def processInputs( lineNamesAr, allCPath, minCov, chrmList, outPre ):
	
	# create this dictionary
	pool = multiprocessing.Pool( processes=NUMPROC )
	print( 'Starting analysis with {:d} processors'.format(NUMPROC) )
	results = [ pool.apply_async( processChrm, args=( lineNamesAr, allCPath, chrm, minCov )) for chrm in chrmList ]
	resultOutput = [ p.get() for p in results ]
	
	print( 'Putting chromosome results together...')
	# put results together
	diffMatrix = []
	posAr = []
	# loop through chromosome outputs
	for r in resultOutput:
		posAr += r[1]
		diffMatrix = addToDifferenceMatrix( diffMatrix, r[0] )
		
	# write output
	outFileSMPStr = outPre + '_cg_smps.tsv'
	outFileDiffMatrix = outPre + '_cg_smps_diff_matrix.tsv'
	print( 'Writing output to {:s} and {:s}...'.format( outFileSMPStr, outFileDiffMatrix ) )
	writeSMPOutput( outFileSMPStr, posAr, minCov, lineNamesAr )
	writeDiffMatrixOutput( outFileDiffMatrix, lineNamesAr, diffMatrix )
	print( 'Done.' )

def parseInputs( argv ):
	
	startInd = 0
	if argv[0].startswith( '-o=' ):
		outPre = argv[0].replace( '-o=', '' )
		startInd += 1
	elif argv[1].startswith( '-o=' ):
		outPre = argv[1].replace( '-o=', '' )
		startInd += 1
	else:
		outPre = 'out'
	
	if argv[0].startswith( '-c=' ):
		chrmList = argv[0].replace( '-c=', '' ).split( ',' )
		startInd += 1
	elif argv[1].startswith( '-c=' ):
		chrmList = argv[1].replace( '-c=', '' ).split( ',' )
		startInd += 1
	else:
		chrmList = CHRMNAMES
	
	try:
		minCov = int( argv[startInd] )
	except ValueError:
		print( 'ERROR: <min_cov> must be an integer' )
		exit()
	
	allCPath = argv[startInd + 1]
	if os.path.isdir( allCPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allCPath ) )
			exit()
	
	lineNamesAr = []
	for i in range(startInd+2, len(argv)):
		lineNamesAr += [ argv[i] ]
	if len( lineNamesAr ) < 2:
		print( 'ERROR: must specify at least 2 lines' )
		exit()
	processInputs( lineNamesAr, allCPath, minCov, chrmList, outPre )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("python3.4 determine_cg_smps.py [-o=outPre][-c=chromosome_names] <min_cov> <allC_path> <line1_name> <line2_name> [lineN_name]*")
	else:
		parseInputs( sys.argv[1:] )
