import sys, math, glob, multiprocessing, subprocess, os

# Usage: heatmap_enrichment_from_bed_allc_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed/allc_file> <input_bed/allc_file>

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95
ADDCOUNT=4

def processInputs( gffFileStr, fpkmFileStr, bedFileStr, inputFileStr, numProc, numBins, numBinsStream, streamSize, percentile ):
	bedName = getSampleName( bedFileStr )
	inputName = getSampleName( inputFileStr )
	outFileStr = 'enrich_{:s}-{:s}_{:d}_{:d}_{:d}.csv'.format( bedName, inputName, numBins, streamSize, int( percentile * 100 ) )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nPercentile for correction: {:.3f}\n\nEnriching    {:s}    over    {:s}\nOutput written to: {:s}\n'.format( numBins, numBinsStream, streamSize, percentile, bedName, inputName, outFileStr ) )
	
	bedName = getSampleName( bedFileStr )
	inputName = getSampleName( inputFileStr )
	outFileStr = 'enrich_{:s}-{:s}_{:d}_{:d}_{:d}.csv'.format( bedName, inputName, numBins, streamSize, int( percentile * 100 ) )
	# read GFF
	print( 'Reading GFF file...' )
	gffDict, chrmFormat = readGFF( gffFileStr )
	# read FPKM
	print( 'Reading FPKM file...' )
	fpkmAr = getFPKM( fpkmFileStr )
	
	# read and process both bed files
	print( 'Processing BED files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(f, gffDict, fpkmAr, numBins, numBinsStream, streamSize, chrmFormat) ) for f in [ bedFileStr, inputFileStr ] ]
	outMats = [ p.get() for p in results ]
	matrixBed = outMats[0]
	matrixInput = outMats[1]
	
	# calculate enrichment
	print( 'Calculating enrichment...' )
	enrichMat, numList = calculateEnrichment( matrixBed, matrixInput )
	# get threshold
	threshold = calculatePercentile( percentile, numList )
	# write output
	print( 'Writing output to {:s}'.format( outFileStr ) )
	writeOutput( outFileStr, enrichMat, threshold )
	print( 'Done' )

def getSampleName( fileStr ):
	# bed file
	if fileStr.endswith( '.bed' ):
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
	else:
		leftIndex = fileStr.rfind('/')
		sampleName = fileStr[leftIndex+1:]
		aIndex = sampleName.find( '_' )
		bIndex = sampleName[aIndex+1:].find( '_' )
		sampleName = sampleName[:bIndex+aIndex+1] + '_mCG'
	return sampleName

def readGFF(gffFileStr):
	
	gffFile = open (gffFileStr, 'r' )
	gffDict = {}
	chrmFormat = None
	
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
		if chrmFormat == None:
			if chrm.isdigit():
				chrmFormat = 'digit'
			elif chrm.startswith( 'Chr0' ):
				chrmFormat = 'zeroed'
			elif chrm.startswith( 'Chr' ):
				chrmFormat = 'chr'
			else:
				chrmFormat = 'asis'
		strand = lineAr[6]
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			# put into geneDict
			gffDict[name] = (chrm, start, end, strand)
	
	gffFile.close()
	return gffDict, chrmFormat

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def getFPKM( fpkmFileStr ):
	'''
		returns the array of genes ordered by fpkm
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
	
	# if exists - check date modified
	if fileExists:
		mTime = os.path.getmtime( fpkmFileStr )
		tTime = os.path.getmtime( fpkmFileStr + ".txt" )
		if tTime < mTime:
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
		if len( lineAr ) == 3:
			# header
			if line.startswith( 'locus' ):
				continue
			fpkm = float( lineAr[2] )
			name = lineAr[4]
			if fpkmDict.get( fpkm ) == None:
				fpkmDict[fpkm] = [n]
			# case 2: fpkm in dict
			else:
				fpkmDict[fpkm] += [name]
		# Cufflinks output
		# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
		# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
		# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
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

def processFile( fileStr, gffDict, fpkmAr, numBins, numBinsStream, streamSize, chrmFormat ):
	isBed = checkBed( fileStr )
	# handle bed files
	if isBed:
		print( '-Reading {:s}...'.format( fileStr ) )
		bedDict, readCounts = readBed( fileStr )
		print( '-Processing {:s}...'.format( fileStr ) )
		outMatrix = processBed( bedDict, gffDict, fpkmAr, numBins, numBinsStream, streamSize, readCounts )
	# handle allc files
	else:
		print( '-Reading {:s}...'.format( fileStr ) )
		allCDict = readAllC( fileStr, chrmFormat )
		print( '-Processing {:s}...'.format( fileStr ) )
		outMatrix = processAllC( allCDict, gffDict, fpkmAr, numBins, numBinsStream, streamSize )
	return outMatrix
		
def checkBed( fileStr ):
	if fileStr.endswith( '.bed' ):
		return True
	else:
		# allC file
		return False

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
	readCount = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		lineAr =line.rstrip().split('\t')
		try:
			curChrm = lineAr[0]
			pos = int( lineAr[1] ) + 1
		except ValueError:
			pass
		else:
			# no dictionary for chrm
			if bedDict.get(curChrm) == None:
				bedDict[curChrm] = {pos:1}
			# dictionary for chrm but position not included
			elif bedDict.get(curChrm).get(pos) == None:
				bedDict[curChrm][pos] = 1
			# dictionary for chrm and position included
			else:
				bedDict[curChrm][pos] += 1
			readCount += 1
	
	bedFile.close()
	return bedDict, readCount

def processBed( bedDict, gffDict, fpkmAr, numBins, numBinsStream, streamSize, readCounts ):
	outMatrix = []
	
	# loop genes in fpkmAr
	for gene in fpkmAr:
		info = gffDict.get(gene)
		if info == None:
			continue
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		outAr = []
		
		# upstream/downstream if needed
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			curUpVarAr = countByRegion( bedDict, chrm, upstream, start, numBinsStream )
			curDownVarAr = countByRegion( bedDict, chrm, end, downstream, numBinsStream )
		else:
			curUpVarAr = []
			curDownVarAr = []
		# gene body part
		curGeneVarAr = countByRegion( bedDict, chrm, start, end, numBins )
		
		# forward strand -> good as is
		if strand == '+':
			outAr = curUpVarAr + curGeneVarAr + curDownVarAr
			
		# reverse strand -> correct
		else:
			# upstream <- reverse downstream
			curDownVarAr.reverse()
			# gene <- reverse gene
			curGeneVarAr.reverse()
			# downstream <- reverse upstream
			curUpVarAr.reverse()
			outAr = curDownVarAr + curGeneVarAr + curUpVarAr
		# add to matrix
		# normalize array for read counts
		outAr2 = [ ( x / float( readCounts ) * 10**9 ) for x in outAr ]
		outMatrix.append( outAr2 )
	return outMatrix

def countByRegion( bedDict, chrm, start, end, numBins ):
	
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print 'binWidth', binWidth
	varAr = [ADDCOUNT] * numBins
	
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
	# correct array for bin width
	outAr = [ x / float(binWidth) for x in varAr ]
	return outAr

def readAllC( fileStr, chrmFormat ):
	'''
		This allC file should be the allC information for all chromosomes in one 
		file not just a single chromosome
	'''
	allCFile = open( fileStr, 'r' )
	allCDict = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if lineAr[3].startswith( 'CG' ):
			# no dictionary for scaffold
			chrm = lineAr[0]
			# handle formatting
			if chrmFormat == 'digit' and chrm.isdigit() == False:
				chrm = chrm.replace( 'Chr', '' )
			elif chrmFormat == 'zeroed' and chrm.isdigit():
				chrm = 'Chr{:02d}'.format( int(chrm) )
			elif chrmFormat == 'chr' and chrm.isdigit():
				chrm = 'Chr{:d}'.format( int( chrm ) )
			if allCDict.get( chrm ) == None:
				allCDict[ chrm ] = {}
			allCDict[chrm][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
	
	allCFile.close()
	return allCDict

def processAllC( allCDict, geneDict, fpkmAr, numBins, numBinsStream, streamSize ):
	outMatrix = []
	#print( allCDict.keys() )
	# loop by gene length
	for gene in fpkmAr:
		info = geneDict.get( gene )
		if info == None:
			continue
		#print( info )
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = methByRegion(allCDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = methByRegion(allCDict, chrm, end, downstream, numBinsStream)
		
		curGeneVarAr = methByRegion(allCDict, chrm, start, end, numBins)
		#print( 'curRepeatVarAr', curRepeatVarAr )
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
		outMatrix.append( outAr )
	return outMatrix

def methByRegion(allCDict, chrm, start, end, numBins):
	''' 
	'''
	# adding one methylated and one unmethylated read per bin
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print ('binWidth', binWidth)
	methAr = [1] * numBins
	totalAr = [2] * numBins
	
	# loop for each bin
	for bin in range(numBins - 1):
		# loop for each position in that bin
		for pos in range(binWidth):
			key = bin * binWidth + pos + start
			#print( key )
			try:
				dictEntry = allCDict.get(chrm).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					methAr[bin]+= dictEntry[0]
					totalAr[bin] += dictEntry[1]
			# chrm not there
			except AttributeError:
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = allCDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				methAr[-1] += dictEntry[0]
				totalAr[-1] += dictEntry[1]
		except AttributeError:
			pass
	z = zip( methAr, totalAr )
	methAr = [ (m,t) for m,t in z ]
	outAr = calculateMethylation( methAr, binWidth )
	return outAr

def calculateMethylation( methAr, binWidth ):
	
	outAr = [-1] * len( methAr )
	for i in range(len(methAr)) :
		if methAr[i][0] == 0 and methAr[i][1] == 0:
			outAr[i] = 0
		elif methAr[i][0] != 0 and methAr[i][1] == 0:
			print( methAr[i] )
		else:
			outAr[i] = (float(methAr[i][0]) / float(methAr[i][1])) / binWidth
	return outAr

def calculateEnrichment( bedMatrix, inputMatrix ):
	
	if len( bedMatrix ) != len( inputMatrix ):
		print( 'ERROR: Matrices are not the same size\nBED Matrix has {:d} genes.\nInput Matrix has {:d} genes.'.format( len(bedMatrix), len(inputMatrix)) )
		exit()
	enrichMatrix = []
	numList = []
	# iterate through genes
	for i in range( len(bedMatrix) ):
		if len(bedMatrix[i]) != len(inputMatrix[i]):
			print( 'ERROR: Matrices are not the same size\n(row {:d})\nBED matrix has {:d} columns.\nInput matrix has {:d} columns.'.format(i, len(bedMatrix[i]),len(inputMatrix[i]) ) )
			exit()
		tmp = [ (bedMatrix[i][j])/(inputMatrix[i][j]) for j in range( len(bedMatrix[i]) ) ]
		enrichMatrix.append( tmp )
		numList += tmp
	return enrichMatrix, numList

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

def writeOutput( outFileStr, outMatrix, threshold ):

	outFile = open( outFileStr, 'w' )
	geneCount = len( outMatrix )
	
	for gene in outMatrix:
		for i in range( 0, len(gene) ):
			if gene[i] > threshold:
				gene[i] = threshold
		outStr = 'G{:08},'.format( geneCount )
		outStrAr = [ '{:.5f}'.format(x) for x in gene ]
		outStr += ','.join( outStrAr )
		outFile.write( outStr + '\n' )
		geneCount -= 1
	outFile.close()

def parseInputs( argv ):
	# Usage: enrichment_from_bed_allc.py [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed/allc_file> <input_bed/allc_file>
	
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	percentile = PERCENTILE
	startInd = 0
	
	for i in range( min(5,len(argv)) ):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
			if numProc > 2:
				print( 'WARNING: no need to use more than 2 processors' )
				numProc = 2
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
		elif argv[i].startswith( '-s=' ):
			try:
				streamSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: upstream/downstream size must be integer' )
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
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStr = argv[startInd+2]
	inputFileStr = argv[startInd+3]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStr, inputFileStr, numProc, numBins, numBinsStream, streamSize, percentile )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: heatmap_enrichment_from_bed_allc_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed/allc_file> <input_bed/allc_file>")
	else:
		parseInputs( sys.argv[1:] )
