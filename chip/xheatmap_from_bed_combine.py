import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: heatmap_from_bed_combine.py [-p=num_proc] [-o=out_id] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed_file> [bed_file]*

NUMPROC=1
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95

def processInputs(gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, outID ):
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nPercentile for correction: {:.3f}\nSamples included: {:s}\n'.format( numBins, numBinsStream, streamSize, percentile, ', '.join(sampleNamesAr) ) )
	print( 'Reading FPKM file...' )
	fpkmAr =  getFPKM( fpkmFileStr )
	print( 'Reading GFF file...' )
	geneDict  = readGFF( gffFileStr )
	
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	results = [ pool.apply_async( processBedFile, args=(bedFileStrAr[i], fpkmAr, geneDict, numBins, numBinsStream, streamSize )) for i in range( len(bedFileStrAr) ) ]
	resMatrix = [ p.get() for p in results ]
	print( 'Combining samples...' )
	combMatrix, combList = combineBed( resMatrix )
	thres = calculatePercentile( percentile, combList )
	print( thres )
	
	outFileStr = 'fpkm_heatmap_' + outID + '.csv'
	info = '#bins_gene:{:d}; bins_stream:{:d}; stream_size:{:d}; percentile:{:.1f}; samples:{:s}'.format( numBins, numBinsStream, streamSize, percentile * 100, ' '.join( sampleNamesAr ) )
	print( 'Writing output to {:s}...'.format( outFileStr) )
	writeOutput( outFileStr, combMatrix, thres, info )
	print( 'Done' )

def getFPKM( fpkmFileStr ):
	'''
		returns the array or genes ordered by fpkm
	'''
	# check if exists
	checkFPKM( fpkmFileStr )
	
	fpkmFile = open( fpkmFileStr + ".txt", 'r' )
	fpkmAr = []
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		fpkmAr += [lineAr[0]]
	fpkmFile.close()
	return fpkmAr

def checkFPKM( fpkmFileStr ):
	'''
		checks of the self-created tab-deliminated fpkm file exists
		it it exists but is outdated, it is recreated
		if it doesn't exist, it is created
	'''
	
	fileExists = os.path.isfile( fpkmFileStr + ".txt" )
	
	# if exists - check date modified and file size
	if fileExists:
		mTime = os.path.getmtime( fpkmFileStr )
		tTime = os.path.getmtime( fpkmFileStr + ".txt" )
		fSize = os.path.getsize( fpkmFileStr + ".txt" )
		if tTime < mTime or fSize < 100:
			print( 'Updating {:s}...'.format(fpkmFileStr + ".txt") )
			createFPKMTextFile( fpkmFileStr )
	# doesn't exist - create
	else:
		createFPKMTextFile( fpkmFileStr )
		print( 'Creating {:s}...'.format(fpkmFileStr + ".txt") )

def createFPKMTextFile( fpkmFileStr ):
	'''
		creates the tab-deliminated FPKM file
		useful so gene order is the same across samples
	'''
	# fpkmDict set up as {fpkm:[gene1,gene2,...], fpkm2:[gene4],...}
	fpkmOutFile = open( fpkmFileStr + ".txt", 'w' )
	fpkmDict = readFPKM( fpkmFileStr )
	
	for key in sorted(fpkmDict.keys(), reverse=True):
		genes = fpkmDict[key]
		# for each gene
		for gene in genes:
			fpkmOutFile.write( "{:s}\t{:.4f}\n".format( gene, key ) )
			
	fpkmOutFile.close()

def readFPKM( fpkmFileStr ):
	'''
		creates a dictionary for genes by fpkm value
		fpkm value is key which points to an array of gene names
		return fpkm dictionary
	'''
	fpkmFile = open( fpkmFileStr, 'r' )
	fpkmDict = {}
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		# Adam's output
		# (0) locus (1) coverage (2) FPKM
		#print( len( lineAr ) )
		if len( lineAr ) == 3:
			# header
			if line.startswith( 'locus' ):
				continue
			fpkm = float( lineAr[2] )
			name = lineAr[0]
			if fpkmDict.get( fpkm ) == None:
				fpkmDict[fpkm] = [name]
			# case 2: fpkm in dict
			else:
				fpkmDict[fpkm] += [name]
		# Cufflinks output
		# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
		# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
		# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
		else:
			if lineAr[12] == 'OK':
				fpkm = float( lineAr[9] )
				name = lineAr[4]
				if name == '-':
					continue
				nameAr = name.split(',')
				for n in nameAr:
					# case 1: fpkm not in fpkmDict
					if fpkmDict.get( fpkm ) == None:
						fpkmDict[fpkm] = [n]
					# case 2: fpkm in dict
					else:
						fpkmDict[fpkm] += [n]
	return fpkmDict
	
def readGFF(gffFileStr):
	
	gffFile = open (gffFileStr, 'r' )
	gffDict = {}
	
	for line in gffFile:
		line = line[:-1]
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		chrm = lineAr[0]
		strand = lineAr[6]
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			# put into geneDict
			gffDict[name] = (chrm, start, end, strand)
	
	gffFile.close()
	return gffDict

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr

def processBedFile( fileStr, fpkmAr, gffDict, numBins, numBinsStream, streamSize ):
	print( 'Reading {:s}...'.format( fileStr ) )
	bedDict, readCounts = readBed( fileStr )
	print( 'Processing {:s}'.format( fileStr ) )
	outMatrix = processBed( bedDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize, readCounts )
	return outMatrix
		
def readBed( bedFileStr ):
	''' takes in the bed file, scaffold we are considering
		returns a dictionary where the position is the key and the value is the 
		number of hits for that position
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	countReads = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		line = line.replace('\n','')
		lineAr =line.split('\t')
		try:
			curChrm = lineAr[0]
			pos = int(lineAr[1]) + 1
		except ValueError:
			pass
		else:
			countReads += 1
			# no dictionary for scaffold
			if bedDict.get(curChrm) == None:
				bedDict[curChrm] = {pos:1}
			# dictionary for scaffold but position not included
			elif bedDict.get(curChrm).get(pos) == None:
				bedDict[curChrm][pos] = 1
			# dictionary for scaffold and position included
			else:
				bedDict[curChrm][pos] += 1
	
	bedFile.close()
	return bedDict, countReads

def processBed( bedDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize, countReads ):
	# fpkmDict organized {fpkm: (gene1, gene2, ...)}
	# geneDict organized {name: (scaffold, start, end, strand)}
	
	outMatrix = []
	# loop by gene fpkm
	for gene in fpkmAr:
		info = geneDict.get( gene )
		# not actually a gene - don't count it
		if info == None:
			continue
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = varByRegion(bedDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = varByRegion(bedDict, chrm, end, downstream, numBinsStream)
		# gene
		curGeneVarAr = varByRegion(bedDict, chrm, start, end, numBins)
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				outAr = curUpVarAr + curGeneVarAr + curDownVarAr
			else:
				outAr = curGeneVarAr
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				outAr = curDownVarAr + curGeneVarAr + curUpVarAr
			else:
				curGeneVarAr.reverse()
				outAr = curGeneVarAr
		# adjust by library size
		outAr = [ float(x)/countReads*1000000 for x in outAr ]
		outMatrix.append( outAr )
	return outMatrix

def varByRegion( bedDict, chrm, start, end, numBins ):
	''' 
		takes in the variant dictionary generated for a gene and the start and 
		end positions of that gene. Note the dictionary is set up so the key is 	
		the position
		on the scaffold and the value is the frequency of that variant
		returns one array with number of variants separated by bin
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print 'binWidth', binWidth
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
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = bedDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				varAr[-1] += dictEntry
		except AttributeError:
			pass
	
	varAr = adjustArrayLength( binWidth, varAr )
	return varAr

def adjustArrayLength( length, inAr ):
	outAr = [0] * len( inAr )
	for i in range( len(inAr) ):
		outAr[i] = float( inAr[i] ) / float( length )
	return outAr

def combineBed( inMatrix ):
	numList = []
	n = len( inMatrix )
	nRow = len( inMatrix[0] )
	nCol = len( inMatrix[0][0] )
	outMatrix = [ [-1] * nCol for x in range( nRow ) ]
	
	# loop through rows
	for i in range( nRow ):
		# loop through through columns
		for j in range( nCol ):
			kSum = 0
			# loop through samples
			for k in range( n ):
				kSum += inMatrix[k][i][j]
			kAve = float(kSum) / n
			outMatrix[i][j] = kAve
			numList += [ kAve ]
		# end for j
	# end for i
	return outMatrix, numList

def calculatePercentile( percentile, numList ):
	
	if percentile == 1:
		return max( numList )
	numList.sort()
	ind = math.ceil( percentile * len( numList ) - 1 )
	try:
		p = numList[ind]
	except IndexError:
		return numList[-1]
	if p == 0:
		print( '***** not percentile corrected *****' )
		return numList[-1]
	return p

def writeOutput( outFileStr, outMatrix, threshold, info ):
	headerAr = ['geneNum'] + [ 'G{:02}'.format(i) for i in range(len(outMatrix[0])) ]

	outFile = open( outFileStr, 'w' )
	outFile.write( info + "\n" + ','.join( headerAr ) + "\n" )
	geneCount = len( outMatrix )
	
	for gene in outMatrix:
		for i in range( 0, len(gene) ):
			if gene[i] > threshold:
				gene[i] = threshold
		outStr = 'G{:08},'.format( geneCount )
		outStrAr = [ '{:.8f}'.format(x) for x in gene ]
		outStr += ','.join( outStrAr )
		outFile.write( outStr + '\n' )
		geneCount -= 1
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	percentile = PERCENTILE
	outID = 'out'
	startInd = 0
	
	for i in range( min(5, len(argv) ) ):
		# number of processors
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
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
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
	# end for
	
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, outID )


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: heatmap_from_bed_combine.py [-p=num_proc] [-o=out_id] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
