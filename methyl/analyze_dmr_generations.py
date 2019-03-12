import sys, math, glob, multiprocessing, subprocess, os, decimal

# Usage: python3.4 analyze_dmr_generations.py [-v] [-p] [-c] [-m=methTypes] [-o=outPre] dmrFile allcPath sample1 sample2 [sampleN]*

MTYPES=['C','CG','CHG','CHH']

def processInputs( dmrFileStr, allCPath, lineNamesAr, mTypes, outPre, steps ):
	# steps: dmr analysis, pvalue analysis, chi-square analysis
	
	# output file names
	valFileStr = outPre + '_dmr_vals.tsv'
	genFileStr = outPre + '_dmr_gen.tsv'
	chiFileStr = outPre + '_dmr_chi.tsv'
	
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	print( 'Methylation Types: {:s}\nSamples included: {:s}\nNumber of DMRs included: {:d}\n'.format( '  '.join(mTypes), '  '.join(lineNamesAr), dmrCount ) )
	print( 'Begin DMR region methylation analysis...' )
	mDict, regionStrAr = computeRegionsMethylation( dmrDict, lineNamesAr, allCPath, mTypes )
	
	if steps[0]:
		print( 'Writing DMR values output to {:s}...'.format( valFileStr ))
		writeValOutput( valFileStr, lineNamesAr, regionStrAr, mDict )
	
	# p-value analysis
	if steps[1]:
		genDict = {}
		# loop through mTypes
		print( 'Begin p-value generation analysis...' )
		for mType in sorted(mDict.keys()):
			print( 'Processing methylation type', mType )
			mAr = generationChange( mDict[mType], regionStrAr, lineNamesAr )
			genDict[mType] = mAr
	
		print( 'Writing generation analysis output to {:s}...'.format( genFileStr ) )
		writeGenOutput( genFileStr, genDict)
	
	# chi-square analysis
	if steps[2]:
		chiDict = {}
		print( 'Beging chi-squared analysis...' )
		for m in sorted(mDict.keys()):
			print( 'Processing methylation type', m  )
			chiDict[m] = processMType( mDict[m] )
		df = len( lineNamesAr ) - 1
		print( 'Writing chi-squared output to {:s}...'.format( chiFileStr ) )
		writeChiOutput( chiFileStr, chiDict, regionStrAr, df)
	print( 'Done' )

def readDMRFile( dmrFileStr ):
	'''
		return dictionary of subset (DMR) regions
		{chr1:[(start,end),(start2,end2)],chr2:[(start,end)],...}
	'''
	dmrFile = open( dmrFileStr, 'r' )
	dmrDict = {}
	dmrCount = 0
	
	for line in dmrFile:
			lineAr = line.rstrip().split()
			chrm = lineAr[0]
			start = int( lineAr[1] )
			end = int( lineAr[2] )
			if dmrDict.get(chrm) == None:
				dmrDict[chrm] = []
			dmrDict[chrm] += [(start, end)]
			dmrCount += 1
	dmrFile.close()
	return dmrDict, dmrCount

def computeRegionsMethylation( dmrDict, lineNamesAr, allCPath, mTypes ):
	
	outDict = {}
	regionStrAr = []
	for m in mTypes:
		outDict[m] = []
	
	# loop through chromosomes
	for chrm in sorted(dmrDict.keys()):
		#						region1		  region2
		# chrmDict[mType] = [ [line1,line2],[line1,line2] ]
		print( 'Processing chromosome {:s}...'.format( chrm ) )
		chrmDict, regionStrs = analyzeChrm( chrm, dmrDict[chrm],lineNamesAr, allCPath, mTypes )
		regionStrAr += regionStrs
		# add to outDict for mTypes
		for m in mTypes:
			outDict[m] = addToMatrix( outDict[m], chrmDict[m] )
	return outDict, regionStrAr

def analyzeChrm( chrm, regionAr, lineNames, allCPath, mTypes ):
	#						region1		  region2
	# outDict[mType] = [ [line1,line2],[line1,line2] ]
	outDict = {}
	regionStrAr = []
	# set up dictionary
	for m in mTypes:
		outDict[m] = [[-1]*len(lineNames) for r in regionAr]
	for r in regionAr:
		rStr = '{:s}:{:d}-{:d}'.format(chrm, r[0], r[1] )
		regionStrAr += [ rStr ]
	# loop through lines
	for j in range(len(lineNames)):
		allCFileStr =os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, lineNames[j], chrm ))
		# ordered by mTypes (rows) then regions (columns)
		lineCAr = analyzeSampleChrm( allCFileStr, regionAr, mTypes )
		# loop through mTypes
		for m in range(len(mTypes)):
			# loop through regions
			for i in range(len(regionAr)):
				#print( 'm',m,'i',i,'j',j)
				#print( outDict )
				#print( lineCAr )
				outDict[mTypes[m]][i][j]= lineCAr[m][i]
	return outDict, regionStrAr

def analyzeSampleChrm( allCFileStr, regionAr, mTypes ):
	
	outAr = [ [-1] * len( regionAr ) for m in mTypes ]
	allCDict = readAllC( allCFileStr, mTypes )
	
	for i in range( len( regionAr ) ):
		rStart = regionAr[i][0]
		rEnd = regionAr[i][1]
		# return array of length mTypes
		cAr = calculateRegionC( rStart, rEnd, allCDict, mTypes )
		for j in range(len(cAr)):
			outAr[j][i] = cAr[j]
	return outAr

def readAllC( allCFileStr, mTypes ):
	allCFile = open( allCFileStr, 'r' )
	print( 'Reading {:s}...'.format( allCFileStr ) )
	allCDict = {}
	
	for m in mTypes:
		allCDict[m] = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		mLineType = findMethylType( lineAr[3] )
		if mLineType in mTypes or lineAr[3] in mTypes:
			allCDict[mLineType][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
		if 'C' in mTypes:
			allCDict['C'][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
	allCFile.close()
	return allCDict

def findMethylType( mc_class ):
	
	if mc_class.startswith( 'CG' ):
		return 'CG'
	elif mc_class.endswith( 'G' ):
		return 'CHG'
	elif mc_class == 'CNN':
		return 'CNN'
	else:
		return 'CHH'

def calculateRegionC( start, end, allCDict, mTypes ):
	
	methC = [0] * len( mTypes )
	totalC = [0] * len( mTypes )
	for pos in range( start, end+1 ):
		for i in range( len(mTypes) ):
			tup = allCDict.get(mTypes[i]).get(pos)
			if tup != None:
				methC[i] += tup[0]
				totalC[i] += tup[1]
	# calculate percentages
	perC = [0] * len(mTypes)
	#print( methC, totalC )
	for j in range( len(mTypes) ):
		perC[j] = (methC[j],totalC[j])
	return perC

def addToMatrix( oldMatrix, newMatrix ):
	# check same number of columns
	if oldMatrix != [] and len(oldMatrix[0]) != len(newMatrix[0] ):
		print( 'ERROR: matrices need to have same number of columns' )
		exit()
	# loop through new matrix
	for region in newMatrix:
		oldMatrix.append( region )
	return oldMatrix

def generationChange( mMatrix, regionStrAr, lineNamesAr ):
	
	# region, line1, line2, change
	outAr = []
	# loop through regions
	for i in range(len(regionStrAr)):
		# loop through generations:
		for j in range(len(lineNamesAr)-1):
			mLine1 = mMatrix[i][j]
			mLine2 = mMatrix[i][j+1]
			mChange = calculateWeightedMethylation(mLine2[0],mLine2[1]) - calculateWeightedMethylation(mLine1[0],mLine1[1])
			mName = lineNamesAr[j].replace('-','.') + '..' +lineNamesAr[j+1].replace('-','.')
			# note: pval is decimal object
			pval = calculatePVal( mLine1, mLine2 )
			try:
				# lpval will also be decimal object
				lpval = decimal.Decimal(-1) * pval.log10()
			except ValueError:
				lpval = -1
			# convert decimals to float
			outAr.append( [ i, regionStrAr[i], mName, mChange, float(pval), float(lpval) ] )
		# for j
	# for i
	return outAr

def calculateWeightedMethylation( countM, countT ):
	if countM == 0 and countT == 0:
		return 0.0
	elif countM != 0 and countT == 0:
		return float( 'inf' )
	else:
		return float( countM ) /float( countT )

def calculatePVal( tup1, tup2 ):
	'''
		tup1 and tup2 are tuples of (num meth, num total) for a given DMR
		this function calculates the p-value as a two-proprortion z-test
		between these two sets of methylation
		return p-value, decimal
	'''
	# set decimal context
	decimal.getcontext().prec = 53
	#print( tup1, tup2)
	# P = [(X0) + (X1)]/[(N0)+(N1)]
	# X is reads methylated, N is total reads
	P = decimal.Decimal( tup1[0] + tup2[0] ) / decimal.Decimal( tup1[1] + tup2[1] )
	p1 = decimal.Decimal(tup1[0] ) / decimal.Decimal(tup1[1])
	p2 = decimal.Decimal(tup2[0] ) / decimal.Decimal(tup2[1])
	#print( P, p1, p2 )
	# z = (p0) - (p1) / sqrt ( (P*(1-P)/(N0))+ (P*(1-P)/(N1)))
	if P == 0 or P == 1:
		return decimal.Decimal(1)
	q = (P*(1-P)/decimal.Decimal(tup1[1])) + (P*(1-P)/decimal.Decimal(tup2[1]))
	z = ( p1 - p2 ) / q.sqrt()
	
	# p(x,u,s) = 1/( s *sqrt(2*pi) ) e^-((x-u)^2/(2*s^2))
	return stdNormDist( abs(z) )
	
def stdNormDist( x ):
	'''
		calculate the p-value of two-proportion z test
		calculates 'cumulative probability density' then adjusts
		for 2-tailed and get value of occurring by chance alone
		
		F(x) = 1/2 ( 1 + erf( x/sqrt(2) ) )
		pvalue = 2 * ( 1 - F(x) )
	'''
	w = x / decimal.Decimal(2.0).sqrt()
	erf = decimal.Decimal(1) - decimal.Decimal( math.erfc(w))
	cdf = ( decimal.Decimal(1) + erf ) / decimal.Decimal(2)
	return decimal.Decimal(2) * ( decimal.Decimal(1) - cdf )

def processMType( mMatrix ):
	
	output = []
	# loop through regions
	for j in range(len(mMatrix)):
		output += [ processRegion( mMatrix[j] ) ]
	return output

def processRegion( inArray ):
	'''
		inArray is array of tuples where each tuple represents the methylated
		and unmethylated reads for a generation
	'''
	# create matrix where rows are methylated counts/unmethylated counts and 
	# columns are lines
	inMatrix = formatMatrix( inArray )
	rowSums, colSums, totalSum = calculateMarginalSums( inMatrix )
	#print( 'observed:' )
	#printChiMatrix(inMatrix, rowSums, colSums, totalSum)
	#print( 'expected:' )
	expectedMatrix = calculatedExpected( rowSums, colSums, totalSum )
	#printChiMatrix(expectedMatrix, rowSums, colSums, totalSum)
	
	chiSquare = computeChiSquared( inMatrix, expectedMatrix )
	return chiSquare

def formatMatrix ( inArray ):
	'''
		convert the matrix with read counts to the chi-square table
		note that inArray contains tuples of (read methylated, total reads)
		so we need to adjust to get reads unmethylated
		row: number reads methylated; number reads unmethylated
		cols: generations included
	'''
	outMat = [ [0] * len( inArray ) for x in range( 2 ) ]
	# loop through tuples
	for i in range(len(inArray)):
		tup = inArray[i]
		outMat[0][i] = tup[0]
		outMat[1][i] = tup[1]-tup[0]
	return outMat

def calculateMarginalSums( inMatrix ):
	'''
		compute row, column, and total marginal sums for a matrix
		used to compute expected values matrix
	'''
	rowSums = [ 0 ] * len(inMatrix)
	colSums = [ 0 ] * len(inMatrix[0])
	# loop through matrix
	for i in range(len(inMatrix)):
		for j in range(len(inMatrix[0])):
			rowSums[i] += inMatrix[i][j]
			colSums[j] += inMatrix[i][j]
	return rowSums, colSums, sum(rowSums)

def calculatedExpected( rowSums, colSums, totalSum ):
	'''
		returns an expected values matrix based on the 
		row, column, and total sums
	'''
	outMat = []
	for i in range(len(rowSums)):
		tmp = [ ( float(rowSums[i]) * x / float(totalSum) ) for x in colSums ]
		outMat.append( tmp )
	return outMat

def computeChiSquared( inMatrix, expectedMatrix ):
	'''
		compute chi-square value given the input matrix and it's
		expected values matrix
	'''
	result = 0
	
	for i in range(len(inMatrix)):
		for j in range(len(inMatrix[0])):
			try:
				result += math.pow( (inMatrix[i][j]-expectedMatrix[i][j]), 2) / float(expectedMatrix[i][j])
			except ZeroDivisionError:
				return 0
	return result

def writeValOutput( outFileStr, lineNamesAr, regionAr, mDict ):
	outFile = open( outFileStr, 'w' )
	# columns: (0) mType (1) position (2) generation (3)
	header = '#mType\tDMR\tregion\tgeneration\tmethylated.reads\ttotal.reads\tweighted.methylation\n'
	outFile.write( header )
	# loop through mTypes
	for mType in sorted( mDict.keys() ):
		mMatrix = mDict[mType]
		# loop through regions
		for i in range(len( regionAr ) ):
			# loop through lines
			for j in range(len(lineNamesAr) ):
				outStr = '{:s}\t{:d}\t{:s}\t{:s}\t{:d}\t{:d}\t{:.6f}\n'.format( mType, i, regionAr[i], lineNamesAr[j], mMatrix[i][j][0],mMatrix[i][j][1], calculateWeightedMethylation( mMatrix[i][j][0], mMatrix[i][j][1]))
				outFile.write( outStr )
			# for j
		# for i
	# for mType
	outFile.close()

def writeGenOutput( outFileStr, outDict):
	# outDict[mType] = [ [regionNum,region,line1,line2,change],...]
	
	outFile = open( outFileStr, 'w' )
	# header
	outFile.write( '#mType\tDMR\tregion\tgeneration\tmeth.change\tpvalue\tlog10.pvalue\n' )
	# loop through mTypes
	for mType in sorted(outDict.keys()):
		# loop through regions
		for region in outDict[mType]:
			#print( region )
			outStr = '{:s}\t{:d}\t{:s}\t{:s}\t{:f}\t{:.10f}\t{:.10f}\n'.format( mType, region[0], region[1], region[2], region[3], region[4], region[5])
			outFile.write( outStr )
	outFile.close()

def writeChiOutput( outFileStr, mDict, regionAr, df):
	
	outFile = open( outFileStr, 'w' )
	# headers 
	outFile.write( '#mType\tDMR\tregion\tchiSquare\tdegrees\n' )
	
	# loop through mType
	for m in sorted(mDict.keys()):
		for j in range(len( regionAr )):
			outStr = '{:s}\t{:d}\t{:s}\t{:.5f}\t{:d}\n'.format(m, j, regionAr[j], mDict[m][j], df )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	'''
	Usage: python3.4 analyze_dmr_generations.py [-v] [-p] [-c] [-o=outPre] [-m=methTypes] dmrFile allcPath sample1 sample2 [sampleN]*
	'''
	
	steps = [0,0,0]
	outPre = 'out'
	startInd = 0
	mTypes = MTYPES
	
	for i in range(min(5,len(argv)-3)):
		if argv[i] == '-v':
			steps[0] = 1
			startInd += 1
		elif argv[i] == '-p':
			steps[1] = 1
			startInd +=1
		elif argv[i] == '-c':
			steps[2] = 1
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			m = argv[i][3:].split(',')
			mTypes = [ x.upper() for x in m ]
			startInd += 1
	# end for
	if sum(steps) == 0:
		steps = [1,1,1]
	
	dmrFileStr = argv[startInd]
	allcPath = argv[startInd+1]
	if os.path.isdir( allcPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allcPath ) )
			exit()
	sampleNamesAr = []
	for i in range(startInd+2, len(argv)):
			sampleNamesAr += [ argv[i] ]
	if len( sampleNamesAr ) < 2:
		print( 'ERROR: must specify at least 2 generations' )
		exit()
	
	processInputs( dmrFileStr, allcPath, sampleNamesAr, mTypes, outPre, steps )

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print ("Usage: python3.4 analyze_dmr_generations.py [-v] [-p] [-c] [-m=methTypes] [-o=outPre] dmrFile allcPath sample1 sample2 [sampleN]*")
		print('-v\toutput dmr methylation values\n-p\toutput p-value analysis\n-c\toutput chi-squared analysis\n[default all on unless one or more is specified]\n-o=outPre\tprefix to use for output files [default:out]\n-m=methTypes\tcomma-separated list of methylation types to include [default:C,CG,CHG,CHH]\ndmfFile\ttab-separated file with DMR coordinates\nallcPath\tpath to the allC files\nsample\tname of samples to be included; order matters for p-value analysis')
	else:
		parseInputs( sys.argv[1:] )
