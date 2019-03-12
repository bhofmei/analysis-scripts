import sys, math, glob, multiprocessing, subprocess, os, datetime

# Usage: metaplot_enrichment_from_bed_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] <gff_file> <bed/allc_file> <input_bed/allc_file>

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000
MINLEN=200
ADDCOUNT=4

# normalization: library size, binWidth, number of features
# enrichment by division

def processInputs( gffFileStr, bedFileStr, inputFileStr, numProc, numBins, numBinsStream, streamSize, typeIncluded, minLen ):

	bedName = getSampleName( bedFileStr )
	inputName = getSampleName( inputFileStr )
	outFileStr = 'meta_enrich_{:s}-{:s}_{:d}_{:d}.csv'.format( bedName, inputName, numBins, streamSize )
	print( '\nNumber of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nFeatures included: {:s}\n\nEnriching    {:s}    over    {:s}\nOutput file: {:s}\n'.format( numBins, numBinsStream, streamSize, ('genes and TEs' if typeIncluded == 'all' else typeIncluded), bedName, inputName, outFileStr ) )
	dt = datetime.datetime.today()
	info = '#created_by:{:s};created:{:d}/{:d}/{:d};num_bins:{:d};num_bins_stream:{:d};stream_size{:d};enriched:{:s};control:{:s}'.format(os.path.basename(__file__), dt.month, dt.day, dt.year, numBins, numBinsStream, streamSize, bedName, inputName )
	
	print( 'Reading GFF file...' )
	geneAr, teAr = readGFF( gffFileStr, typeIncluded, minLen )
	if geneAr != None:
		info += ';num_genes:{:d}'.format( len(geneAr) )
		print( '\tnumber of genes: {:d}'.format( len(geneAr)) )
	if teAr != None:
		info += ';num_tes:{:d}'.format( len(teAr))
		print( '\tnumber of repeats/TEs: {:d}'.format( len(teAr)) )
	
	pool = multiprocessing.Pool( processes=numProc )
	fileAr = [ bedFileStr, inputFileStr ]
	#print( fileAr )
	
	results = [ pool.apply_async(readBed, args=(f, True) ) for f in fileAr ]
	bedDicts = [ p.get() for p in results ]
	bedDict = bedDicts[0]
	inputDict = bedDicts[1]
	
	featMat = [ geneAr, teAr ]
	featStr = [ 'genes', 'TEs/repeats' ]
	
	results = [ pool.apply_async(processFeature, args=(featMat[i], bedDict, inputDict, numBins, numBinsStream, streamSize, featStr[i] )) for i in range(len(featStr)) ]
	metaMatrix = [ p.get() for p in results ]
	geneMetaAr = metaMatrix[0]
	teMetaAr = metaMatrix[1]
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, geneMetaAr, teMetaAr, info )

def getSampleName( fileStr ):
	# bed file
	leftIndex = fileStr.rfind('/')
	rightIndex = fileStr.rfind('.')
	sampleName = fileStr[leftIndex+1:rightIndex]
	return sampleName

def readGFF( gffFileStr, typeIncluded, minLen ):
	# typeIncluded = 'all' | 'genes' | 'tes'
	if typeIncluded == 'all':
		geneAr = []
		teAr = []
	elif typeIncluded == 'genes':
		geneAr = []
		teAr = None
	elif typeIncluded == 'tes':
		geneAr = None
		teAr = []
	
	gffFile = open (gffFileStr, 'r' )
	for line in gffFile:
		line = line[:-1]
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		if ( end - start + 1 ) <= minLen:
			continue
		chrm = lineAr[0]
		strand = lineAr[6]
		
		if typeIncluded in [ 'all', 'tes' ] and lineAr[2] in ['transposable_element', 'transposable_element_gene', 'repeat']:
			teAr += [ (chrm, start, end, strand) ]
		if typeIncluded in [ 'all', 'genes' ] and lineAr[2] == 'gene':
			geneAr+= [(chrm, start, end, strand)]
	# end for line
	gffFile.close()
	return geneAr, teAr
	
def processFeature( featAr, bedDict, inputDict, numBins, numBinsStream, streamSize, strID ):
	
	if featAr == None:
		return None
	
	print( 'Processing {:s}...'.format( strID ) )
	# process bedDict
	print( '-Creating value matrices for', strID )
	bedMatrix = processBed( bedDict, featAr, numBins, numBinsStream, streamSize, bedDict['count'] )
	#print( bedMatrix[0] )
	inputMatrix = processBed( inputDict, featAr, numBins, numBinsStream, streamSize, inputDict['count'] )
	#print( inputMatrix[0] )
	
	print( '-Computing enrichment for', strID )
	# calculate enrichment
	enrichMat = calculateEnrichment( bedMatrix, inputMatrix )
	#print( enrichMat[0] )
	
	print( '-Combining enrichment values for', strID )
	# sum for feature - adjust by number of features
	metaAr = calculateMeta( enrichMat, len(featAr ) )
	#print( metaAr )
	return metaAr
	
def readBed( bedFileStr, unused ):
	''' takes in the bed file, scaffold we are considering
		returns a dictionary where the position is the key and the value is the 
		number of hits for that position
	'''
	print( 'Reading {:s}...'.format( bedFileStr ) )
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
	bedDict['count'] = countReads
	return bedDict

def processBed( bedDict, featAr, numBins, numBinsStream, streamSize, readCounts ):
	outMatrix = []
	for feat in featAr:
		chrm = feat[0]
		start = feat[1]
		end = feat[2]
		strand = feat[3]
		outAr = []
		
		# upstream/downstream
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
		# normalize array for read counts
		outAr2 = [ ( float(x) / float( readCounts ) * 10**9 ) for x in outAr ]
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

def calculateEnrichment( bedMatrix, inputMatrix ):
	if len( bedMatrix ) != len( inputMatrix ):
		print( 'ERROR: Matrices are not the same size\nBED Matrix has {:d} genes.\nInput Matrix has {:d} genes.'.format( len(bedMatrix), len(inputMatrix)) )
		exit()
	enrichMatrix = []
	# iterate through genes
	for i in range( len(bedMatrix) ):
		if len(bedMatrix[i]) != len(inputMatrix[i]):
			print( 'ERROR: Matrices are not the same size\n(row {:d})\nBED matrix has {:d} columns.\nInput matrix has {:d} columns.'.format(i, len(bedMatrix[i]),len(inputMatrix[i]) ) )
			exit()
		tmp = [ (bedMatrix[i][j])/(inputMatrix[i][j]) for j in range( len(bedMatrix[i]) ) ]
		enrichMatrix.append( tmp )
	return enrichMatrix

def calculateMeta( enrichMat, numFeat ):
	
	outAr = [0] * len(enrichMat[0])
	
	# loop through features
	for i in range(len(enrichMat)):
		for j in range(len(enrichMat[0])):
			outAr[j] += enrichMat[i][j]
	# adjust by num features
	adjustAr = [ float(x) / numFeat for x in outAr ]
	return adjustAr

def writeOutput( outFileStr, geneMetaAr, teMetaAr, info ):
	
	outFile = open( outFileStr, 'w' )
	header = '#feature,bin,value\n{:s}\n'.format( info )
	outFile.write( header )
	
	if geneMetaAr != None:
		for i in range(len(geneMetaAr)):
			outStr = 'genes,{:d},{:f}\n'.format( i+1, geneMetaAr[i] )
			outFile.write( outStr )
	if teMetaAr != None:
		for i in range(len(teMetaAr)):
			outStr = 'tees,{:d},{:f}\n'.format( i+1, teMetaAr[i] )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	# Usage: metaplot_enrichment_from_bed_allc.py [-a|-g|-te] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed/allc_file> <input_bed/allc_file>
	
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	minLen = MINLEN
	typeIncluded = None
	startInd = 0
	
	for i in range( min(9,len(argv)-3) ):
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
		elif argv[i].startswith( '-m=' ):
			try:
				minLen = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum repeat length must be integer' )
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
		elif argv[i].startswith( '-s=' ):
			try:
				streamSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: upstream/downstream size must be integer' )
				exit()
		elif argv[i] == '-a':
			if typeIncluded != None:
				print( 'ERROR: cannot specify -a with {:s}'.format( ('-g' if typeIncluded == 'genes' else '-te' ) ) )
				exit()
			typeIncluded = 'all'
			startInd += 1
		elif argv[i] == '-g':
			if typeIncluded != None:
				print( 'ERROR: cannot specify -g with {:s}'.format( ('-a' if typeIncluded == 'all' else '-te' ) ) )
				exit()
			typeIncluded = 'genes'
			startInd += 1
		elif arg[i] == '-te':
			if typeIncluded != None:
				print( 'ERROR: cannot specify -te with {:s}'.format( ('-a' if typeIncluded == 'all' else '-g' ) ) )
				exit()
			typeIncluded = 'tes'
			startInd += 1
	# end for
	if typeIncluded == None:
		typeIncluded = 'all'

	gffFileStr = argv[startInd]
	bedFileStr = argv[startInd+1]
	inputFileStr = argv[startInd+2]
	
	processInputs( gffFileStr, bedFileStr, inputFileStr, numProc, numBins, numBinsStream, streamSize, typeIncluded, minLen )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: metaplot_enrichment_from_bed_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] <gff_file> <bed/allc_file> <input_bed/allc_file>")
	else:
		parseInputs( sys.argv[1:] )
