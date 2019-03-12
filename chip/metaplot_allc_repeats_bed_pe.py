import sys, math, glob, multiprocessing, subprocess, os, datetime

# Usage: python3.4 metaplot_allc_repeats_bed_pe.py [-m=meth_types] [-t=meth_thresholds] [-b=num_bins[,num_bins_stream]] [-s=stream_size]  [-o=out_id] [-p=num_proc] [-r=min_repeat_len] <repeat_gff> <allc_file> <bed_file> [bed_file]*


# Steps
# 1. Repeat array (size limit by bins)
# 2. allC dictionary (methylation type limit)
# 3. Repeat-allc dictionary (limit by methylation threshold)
# 4. bed dictionary
# 5. Compute bed values for each repeat-allc set
# Note: normalization by: library size, repeat length/bin width, number of repeats per methylation type
# 6. Output

MTYPES=['C','CG','CHG','CHH']
MTHRESH=[0.1,0.3,0.2,0.1]
NUMBINS=20
NUMBINSSTREAM=4
STREAMSIZE=200
NUMPROC=1
MINLEN=200

def processInputs( repeatFileStr, allcFileStr, bedFileStrAr, mTypes, mThresholds, numBins, numBinsStream, streamSize, outID, numProc, minLen ):
	
	# check outID
	if outID == None and len( bedFileStrAr ) == 1:
		bName = os.path.basename( bedFileStrAr[0] )
		rInd = bName.rfind( '.' )
		outID = bName[:rInd]
	elif outID == None:
		outID = 'out'
	
	sampleNamesAr = getSampleNames( bedFileStrAr )
	
	# read repeat file
	print( 'Reading {:s}...'.format( repeatFileStr ) )
	repeatAr = readRepeatGFF( repeatFileStr, numBins )
	# read allc file
	print( 'Reading {:s}...'.format( allcFileStr ) )
	allcDict = readAllC( allcFileStr, mTypes )
	# trim repeat regions by methylation level
	print( 'Filtering repeat regions based on methylation thresholds...' )
	repeatDict = filterRepeats( repeatAr, allcDict, mTypes, mThresholds )
	
	mCounts = []
	for m in mTypes:
		mCounts += [ str( len(repeatDict[m]) ) ]
		
	dt = datetime.datetime.today()
	mThresholdStr = [ '{:.2f}'.format( i ) for i in mThresholds ]
	info = '#created_by:{:s};created:{:d}/{:d}/{:d};num_bins:{:d};num_bins_stream:{:d};stream_size:{:d};meth_types:{:s};meth_thresholds:{:s};meth_feat_counts:{:s}'.format( os.path.basename(__file__), dt.month, dt.day, dt.year, numBins, numBinsStream, streamSize, ','.join(mTypes), ','.join( mThresholdStr ), ','.join( mCounts ) )
	
	# process BED for each mType
	print( 'Begin processing BED files with {:d} processors...'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(bedFileStr, repeatDict, numBins, numBinsStream, streamSize, mTypes) ) for bedFileStr in bedFileStrAr ]
	outSuperMatrix = [ p.get() for p in results ]
	# outSuperMatrix [bedFile][mType][bins]
	
	outFileStr = 'meta_meth_repeats_{:s}.tsv'.format( outID )
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, outSuperMatrix, info, mTypes, sampleNamesAr )
	print( 'Done.' )

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr

def readRepeatGFF( repeatFileStr, numBins ):
	
	repeatAr = []
	repeatFile = open( repeatFileStr, 'r' )
	
	for line in repeatFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) source (2) feature (3) start (4) end (5) ?
		# (6) strand (7) ? (8) notes
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		fLen = end - start + 1
		if fLen < numBins:
			continue
		chrm = lineAr[0]
		strand = lineAr[6] 
		repeatAr += [ (chrm, start, end, strand) ]
	repeatFile.close()
	return repeatAr

def readAllC( allCFileStr, mTypes ):
	allCFile = open( allCFileStr, 'r' )
	allCDict = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		chrm = lineAr[0]
		if allCDict.get( chrm ) == None:
			allCDict[chrm] = {}
		pos = int( lineAr[1] )
		mLineType = findMethylType( lineAr[3] )
		if mLineType in mTypes or lineAr[3] in mTypes:
			allCDict[chrm][pos] = ( mLineType, int(lineAr[4]), int( lineAr[5]) )
		#if 'C' in mTypes:
			#allCDict[chrm][pos] = ( 'C', int(lineAr[4]), int( lineAr[5]) )
	allCFile.close()
	return allCDict

def findMethylType( mc_class ):
	
	if mc_class.startswith( 'CG' ):
		return 'CG'
	elif mc_class.startswith( 'C' ) and mc_class.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'

def filterRepeats( repeatAr, allcDict, mTypes, mThresholds ):
	
	repeatDict = {}
	for m in mTypes:
		repeatDict[m] = []
	
	# loop through repeatAr
	for repeat in repeatAr:
		chrm = repeat[0]
		start = repeat[1]
		end = repeat[2]
		cDict = allcDict.get(chrm)
		if cDict == None:
			continue
		wMeth = calculateRegionMethylation( cDict, start, end, mTypes )
		# check vs thresholds
		for i in range(len(wMeth)):
			if wMeth[i] >= mThresholds[i]:
				repeatDict[ mTypes[i] ] += [ repeat ]
	# end for repeat
	return repeatDict

def calculateRegionMethylation( allCDict, start, end, mTypes ):
	# allCDict is specified by chromosome already
	
	methCount = [0] * len( mTypes )
	totalCount = [0] * len( mTypes )
	for pos in range( start, end + 1 ):
		dictEntry = allCDict.get( pos )
		if dictEntry != None:
			if dictEntry[0] in mTypes:
				mInd = mTypes.index( dictEntry[0] )
				methCount[mInd] += dictEntry[1]
				totalCount[mInd] += dictEntry[2]
			if 'C' in mTypes:
				mInd = mTypes.index( 'C' )
				methCount[mInd] += dictEntry[1]
				totalCount[mInd] += dictEntry[2]
		# end if dictEntry != None
	# end for
	outAr = [ calculateWeightedMethylation( methCount[i], totalCount[i] ) for i in range(len(mTypes)) ]
	return outAr

def calculateWeightedMethylation( mCount, tCount ):
	
	if mCount == 0 and tCount == 0:
		return 0.0
	elif mCount == 0 and tCount != 0:
		return float( 'nan' )
	else:
		return float( mCount ) / float( tCount )

def processFile( bedFileStr, repeatDict, numBins, numBinsStream, streamSize, mTypes ):
	print( 'Reading {:s}...'.format( bedFileStr ) )
	bedDict, readCounts = readBed( bedFileStr )
	outMatrix = []
	
	for m in mTypes:
		print( 'Processing methylation type {:s} for {:s}...'.format( m, bedFileStr ) )
		tmpAr = processBed( repeatDict[m], bedDict, readCounts, numBins, numBinsStream, streamSize )
		outMatrix.append( tmpAr )
	print( 'Finished processing {:s}.'.format( bedFileStr ) )
	del( bedDict )
	return outMatrix
		
def readBed( bedFileStr ):
	'''
		creates a dictionary with each scaffold and those dictionaries are a 
		dictionary
		for the positions of coordinates with frequency count from the bed file
		{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
		return the dictionary
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	readCounts = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		lineAr =line.rstrip().split('\t')
		try:
			curScaf = lineAr[0]
			pos = int( lineAr[1] ) + 1
		except ValueError:
			pass
		else:
			# no dictionary for scaffold
			if bedDict.get(curScaf) == None:
				bedDict[curScaf] = {pos:1}
			# dictionary for scaffold but position not included
			elif bedDict.get(curScaf).get(pos) == None:
				bedDict[curScaf][pos] = 1
			# dictionary for scaffold and position included
			else:
				bedDict[curScaf][pos] += 1
			readCounts += 1
	
	bedFile.close()
	return bedDict, readCounts

def processBed( repeatAr, bedDict, countReads, numBins, numBinsStream, streamSize ):
	# repeatAr organized [(scaffold, start, end, strand), ...]
	
	totalRepeatAr = [0] * numBins
	totalUpAr = []
	totalDownAr = []
	
	if streamSize != 0:
		totalUpAr = [0] * numBinsStream
		totalDownAr = [0] * numBinsStream
	repeatCount = len( repeatAr )
	
	# loop by gene length
	for repeat in repeatAr:
		chrm = repeat[0]
		start = repeat[1]
		end = repeat[2]
		strand = repeat[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = varByRegion(bedDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = varByRegion(bedDict, chrm, end, downstream, numBinsStream)
		
		# repeat region
		curRepeatVarAr = varByRegion(bedDict, chrm, start, end, numBins)
		#print( curGeneVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalRepeatAr = addToArray( totalRepeatAr, curRepeatVarAr )
			else:
				totalRepeatAr = addToArray( totalRepeatAr, curRepeatVarAr )
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curRepeatVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				#outAr = curDownVarAr + curRepeatVarAr + curUpVarAr
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalRepeatAr = addToArray( totalRepeatAr, curRepeatVarAr )
			else:
				curRepeatVarAr.reverse()
				totalRepeatAr = addToArray( totalRepeatAr, curRepeatVarAr )
		# end if-else strand
	# end for gene
	outAr = totalUpAr + totalRepeatAr + totalDownAr
	#print( outAr )
	outAr = [ float(x)/( countReads * repeatCount ) * 1000000000 for x in outAr ]
	return outAr

def varByRegion(bedDict, chrm, start, end, numBins):
	''' 
		takes in the variant dictionary generated for a gene and the start and 
		end positions of that gene. 
		Note the dictionary is set up so the key is the position
		on the scaffold and the value is the frequency of that variant
		returns one array with number of variants separated by bin
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	if binWidth < 1:
		return False
	varAr = [0] * numBins
	
	# loop for each bin
	for bin in range(numBins - 1):
		# loop for each position in that bin
		for pos in range(binWidth):
			key = bin * binWidth + pos + start
			try:
				dictEntry = bedDict.get(chrm).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					varAr[bin] += dictEntry
			except AttributeError:
				#print( 'ERROR with', chrm )
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = bedDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				varAr[-1] += dictEntry
		except AttributeError:
			#print( 'ERROR with', chrm )
			pass
	varAr = adjustArrayLength( binWidth, varAr )
	return varAr

def adjustArrayLength( length, inAr ):
	outAr = [0] * len( inAr )
	for i in range( len(inAr) ):
		outAr[i] = float( inAr[i] ) / float( length )
	return outAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	return oldAr

def writeOutput( outFileStr, outSuperMatrix, info, mTypes, sampleNamesAr ):
	outFile = open( outFileStr, 'w' )
	header = 'sample\tmType\tbin\tvalue'
	outFile.write( '{:s}\n{:s}\n'.format( header, info ) )
	# outSuperMatrix [bedFile][mType][bins]
	# loop through bedFile
	for i in range(len(sampleNamesAr)):
		# loop through mTypes
		for j in range(len(mTypes)):
			# loop through bins
			outStrAr = [ '{:s}\t{:s}\t{:d}\t{:f}'.format( sampleNamesAr[i], mTypes[j], k+1, outSuperMatrix[i][j][k] ) for k in range(len(outSuperMatrix[i][j] )) ]
			outFile.write( '{:s}\n'.format( '\n'.join( outStrAr ) ) )
	# end for i
	outFile.close()

def parseInputs( argv ):
	
	mTypes = MTYPES
	mThresholds = None
	numBins = NUMBINS
	numBinsStream = NUMBINSSTREAM
	streamSize = STREAMSIZE
	outID = None
	numProc = NUMPROC
	minLen = MINLEN
	startInd = 0
	
	for i in range(min(7,len(argv)-3)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-r=' ):
			try:
				minLen = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum repeat length must be integer' )
				exit()
		elif argv[i].startswith( '-m=' ):
			mTypes = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			tmpT = argv[i][3:].split(',')
			try:
				tmpM = [ float(x) for x in tmpT ]
				mThresholds = [ ( x if x <= 1 else x / 100 ) for x in tmpM ]
				startInd += 1
			except ValueError:
				print( 'ERROR: check that methylation thresholds are numeric and comma-separated' )
				exit()
		elif argv[i].startswith( '-s=' ):
			try:
				streamSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: upstream/downstream size must be integer' )
				exit()
		elif argv[i].startswith( '-b=' ):
			try:
				numBinsAr = argv[i][3:].split(',')
				numBins = int( numBinsAr[0] )
				if len( numBinsAr ) == 2:
					numBinsStream = int( numBinsAr[1] )
				else:
					numBinsStream = int( numBinsAr[0] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of bins must be integer' )
				exit()
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
	# end for
	
	if mThresholds == None:
		tmpAr = []
		for i in range(len(mTypes)):
			try:
				mInd = MTYPES.index( mTypes[i] )
				tmpAr += [ MTHRESH[mInd] ]
			except ValueError:
				print( 'ERROR: custom or incorrect methylation type specified...specify methylation thresholds' )
				exit()
		mThresholds = tmpAr
	elif len(mTypes) != len( mThresholds):
		print( "ERROR: number of methylation thresholds doesn't match number of methylation types. Specify the methylation types and specify the same number of thresholds." )
		exit()
	repeatFileStr = argv[startInd]
	allcFileStr = argv[startInd + 1]
	bedFileStrAr = []
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( repeatFileStr, allcFileStr, bedFileStrAr, mTypes, mThresholds, numBins, numBinsStream, streamSize, outID, numProc, minLen )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 metaplot_allc_repeats_bed_pe.py [-m=meth_types] [-t=meth_thresholds] [-b=num_bins[,num_bins_stream]] [-s=stream_size]  [-o=out_id] [-p=num_proc] <repeat_gff> <allc_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
