import sys, math, glob, multiprocessing, subprocess, os

'''
	reports the weighted methylation over a set of input regions
	for the input samples specified
'''

# Usage: python3.4 dmr_methylation.py [-r] [-g=genome_methylation] [-p=num_proc] [-m=meth_type] [-o=outPre] <dmr_file> <allC_path> <sample_name> [sample_name]*

NUMPROC=2
MTYPES=['C','CG','CHG','CHH']

def processInputs( dmrFileStr, allCPath, lineNamesAr, mTypes, numProc, outPre, readCounts, genomeFile ):
	
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	print( 'Methylation types: {:s}\nSamples included: {:s}\nNumber of DMRs: {:d}\nGenome methylation correction: {:s}\n'.format( '  '.join(mTypes), '  '.join(lineNamesAr),dmrCount, str(genomeFile != None) ) )
	
	if genomeFile != None:
		correctDict = readGenomeFile( genomeFile, lineNamesAr, mTypes )
		print( 'Read genome correction file.' )
	else:
		correctDict = None
	
	mDict, regionStrAr = processDMRs( dmrDict, lineNamesAr, allCPath, mTypes, readCounts, correctDict )
	print( 'Methylation of regions for all samples calculated' )
	outFileStr = outPre + '_parsed_dmrs.csv'
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, lineNamesAr, regionStrAr, mDict, readCounts, genomeFile )
	print( 'Done.')

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

def readGenomeFile( genomeFileStr, sampleNamesAr, mTypes ):
	
	genomeFile = open( genomeFileStr, 'r' )
	correctionDict = {}	# key is sampleName, value is dict; key is mType, value is weighted methylation
	isFirst = True
	mTypeInd = [-1]*len(mTypes)
	
	for line in genomeFile:
		lineAr = line.rstrip().split('\t')
		# (0) sample (1)+ mTypes
		# header
		if isFirst:
			isFirst = False
			for i in range(len(mTypes)):
				try:
					mTypeInd[i] = lineAr.index(mTypes[i])
				except ValueError:
					print( 'ERROR: genome file does not contain all of the methylation types' )
					exit()
		# other lines
		else:
			samp = lineAr[0]
			correctionDict[samp] = {}
			for j in range(len(mTypes)):
				correctionDict[samp][mTypes[j]] = float( lineAr[ mTypeInd[j] ] )
	genomeFile.close()
	# check all samples included
	for name in sampleNamesAr:
		if correctionDict.get(name) == None:
			print( 'ERROR: genome file does not contain all the samples' )
			exit()
	return correctionDict

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
	elif mc_class.startswith( 'C' ) and mc_class.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'

def calculateRegionC( start, end, allCDict, readCounts, mTypes, correctMDict ):
	
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
	for j in range( len(mTypes) ):
		if readCounts:
			perC[j] = '{:d},{:d}'.format(methC[j],totalC[j])
		elif methC[j] == 0 and totalC[j] == 0 and correctMDict == None:
			perC[j] = '0.0'
		elif methC[j] == 0 and totalC[j] == 0:
			perC[j] = '0.0,0.0'
		elif methC[j] != 0 and totalC == [0] and correctMDict == None:
			perC[j] = '-1'
		elif methC[j] != 0 and totalC == [0]:
			perC[j] = '-1,-1'
		elif correctMDict == None:
			perC[j] = str(float(methC[j]) / float(totalC[j]))
		else:
			perC[j] = str(float(methC[j]) / float(totalC[j])) + ',' +str( float(methC[j]) / float(totalC[j]) / correctMDict[mTypes[j]])
	return perC

def analyzeChrm( chrm, regionAr, lineNames, allCPath, mTypes, readCounts, correctDict ):
	#						region1		  region2
	# outDict[mType] = [ [line1,line2],[line1,line2] ]
	print( 'Analyzing chromosome {:s}...'.format( chrm ) )
	outDict = {}
	regionStrAr = []
	# set up dictionary
	for m in mTypes:
		outDict[m] = [[-1]*len(lineNames) for r in regionAr]
	# loop through lines
	for j in range(len(lineNames)):
		if lineNames[j].find(',') == -1:
			allCFileStr =os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, lineNames[j], chrm ) )
			# ordered by mTypes (rows) then regions (columns)
			if correctDict == None:
				p = None
			else:
				p = correctDict[lineNames[j]]
			lineCAr = analyzeSampleChrm( allCFileStr, regionAr, mTypes, readCounts, p )
			# loop through mTypes
			for m in range(len(mTypes)):
				# loop through regions
				for i in range(len(regionAr)):
					#print( 'm',m,'i',i,'j',j)
					#print( outDict )
					#print( lineCAr )
					outDict[mTypes[m]][i][j]= lineCAr[m][i]
		# need to average samples together
		else:
			lNamesAr = lineNames[j].split(',')
			linesCs = []
			for lName in lNamesAr:
				allCFileStr = os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, lName, chrm ) )
				if correctDict == None:
					p = None
				else:
					p = correctDict[lineNames[j]]
				tmpAr = analyzeSampleChrm( allCFileStr, regionAr, mTypes, readCounts, p )
				linesCs.append( tmpAr )
			# loop through mTypes
			for m in range(len(mTypes)):
				# loop through regions
				for i in range(len(regionAr)):
					if readCounts or correctDict == None:
						print( 'ERROR: cannot use read count option or genome correction with joined lines' )
						exit()
					s = 0
					for k in range(len(lNamesAr)):
						s += float(linesCs[k][m][i])
					outDict[mTypes[m]][i][j] = str( s / len(lNamesAr) )
				# end i
			# end m
		# end else
	# end j
						
	print( 'Finished analyzing chromosome {:s}'.format( chrm ) )
	return outDict

def analyzeSampleChrm( allCFileStr, regionAr, mTypes, readCounts, correctMDict ):
	
	outAr = [ [-1] * len( regionAr ) for m in mTypes ]
	allCDict = readAllC( allCFileStr, mTypes )
	
	for i in range( len( regionAr ) ):
		rStart = regionAr[i][0]
		rEnd = regionAr[i][1]
		# return array of length mTypes
		cAr = calculateRegionC( rStart, rEnd, allCDict, readCounts, mTypes, correctMDict )
		for j in range(len(cAr)):
			outAr[j][i] = cAr[j]
	return outAr

def processDMRs( dmrDict, lineNamesAr, allCPath, mTypes, readCounts, correctDict ):
	
	outDict = {}
	for m in mTypes:
		outDict[m] = []
	regionStrAr = []
	dictKeys = sorted( dmrDict.keys() )
	# fill regionStrAr
	for chrm in dictKeys:
		regionAr = dmrDict[chrm]
		for region in regionAr:
			regionStrAr += [ '{:s}:{:d}-{:d}'.format( chrm, region[0],region[1] ) ]
	# loop through chromosomes
	pool = multiprocessing.Pool( processes=NUMPROC )
	print( 'Begin processing with {:d} processors'.format( NUMPROC ) )
	results = [ pool.apply_async( analyzeChrm, args=(chrm, dmrDict[chrm], lineNamesAr, allCPath, mTypes, readCounts, correctDict) ) for chrm in dictKeys ]
	rDicts = [ p.get() for p in results ]
	
	# add to outDict
	for i in range(len(rDicts)):
		for m in mTypes:
			outDict[m] = addToMatrix( outDict[m], rDicts[i][m] )
	return outDict, regionStrAr

def writeOutput( outFileStr, lineNamesAr, regionAr, mDict, readCounts, genomeFile ):
	outFile = open( outFileStr, 'w' )
	# columns: (0) mType (1) position (2) generation (3)
	if readCounts:
		header = '#mType,DMR,region,generation,methylated_reads,total_reads\n'
	elif genomeFile != None:
		header = '#mType,DMR,region,generation,weighted_methylation,correct_methylation\n'
	else:
		header = '#mType,DMR,region,generation,weighted_methylation\n'
	outFile.write( header )
	# loop through mTypes
	for mType in sorted( mDict.keys() ):
		mMatrix = mDict[mType]
		# loop through regions
		for i in range(len( regionAr ) ):
			# loop through lines
			for j in range(len(lineNamesAr) ):
				outStr = '{:s},{:d},{:s},{:s},{:s}\n'.format( mType, i, regionAr[i], lineNamesAr[j].replace(',',';'), mMatrix[i][j] )
				outFile.write( outStr )
			# for j
		# for i
	# for mType
	outFile.close()

def addToMatrix( oldMatrix, newMatrix ):
	# check same number of columns
	if oldMatrix != [] and len(oldMatrix[0]) != len(newMatrix[0] ):
		print( 'ERROR: matrices need to have same number of columns' )
		exit()
	# loop through new matrix
	for region in newMatrix:
		oldMatrix.append( region )
	return oldMatrix

def parseInputs( argv ):
	# Usage: python3.4 dmr_methylation.py [-r] [-g=genome_methylation] [-p=num_proc] [-m=meth_type] [-o=outPre] <dmr_file> <allC_path> <sample_name> [sample_name]*

	mTypes = MTYPES
	numProc = NUMPROC
	outPre = 'dmr_out'
	readCounts = False
	genomeFile = None
	startInd = 0
	
	for i in range( min(6,len(argv)-3) ):
		if argv[i] == '-r':
			readCounts = True
			startInd += 1
		elif argv[i].startswith('-m'):
			mTypes = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-g=' ):
			genomeFile = argv[i][3:]
			startInd += 1
	# end for
	dmrFileStr = argv[startInd]
	
	allCPath = argv[startInd+1]
	if os.path.isdir( allCPath ) == False:
		print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allCPath ) )
		exit()
		
	lineNamesAr = []
	for i in range(startInd+2, len(argv)):
			lineNamesAr += [ argv[i] ]
	processInputs( dmrFileStr, allCPath, lineNamesAr, mTypes, numProc, outPre, readCounts, genomeFile )
	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 dmr_methylation.py [-r] [-f=genome_methylation] [-m=meth_type] [-p=num_proc] [-o=outPre] <dmr_file> <allC_path> <sample_name> [sample_name,sample_name] [sample_name]*")
	else:
		parseInputs( sys.argv[1:] )
