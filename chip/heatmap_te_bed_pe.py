import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 heatmap_te_bed_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <bed_file> [bed_file]*

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95

def processInputs( gffFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile ):
	
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nPercentile for correction: {:.3f}\nSamples included: {:s}\n'.format( numBins, numBinsStream, streamSize, percentile, ', '.join(sampleNamesAr) ) )
	print( 'Reading GFF file...' )
	gffDict = readGFF( gffFileStr )
	gffAr = processGFFDict( gffDict )
	info = "#from_script:heatmap_te_bed_pe.py; num_bins:{:d}; num_bins_stream:{:d}; stream_size:{:d}; percentile:{:.1f}; gff_file:{:s}".format( numBins, numBinsStream, streamSize, percentile, os.path.basename( gffFileStr ) )
	
	outFileStrAr = [ '{:s}_te_{:d}_{:d}_{:d}.csv'.format( sampleName, streamSize, numBins, int(percentile*100) ) for sampleName in sampleNamesAr ]
	
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	results = [ pool.apply_async( processFile, args=(bedFileStrAr[i], outFileStrAr[i], gffAr, numBins, numBinsStream, streamSize, percentile, info )) for i in range( len(bedFileStrAr) ) ]
	fin = [ p.get() for p in results ]
	if sum(fin) != len( bedFileStrAr ):
		print( 'ERROR: not enough files written' )
		exit()
	print( 'Done' )


def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr

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
		if lineAr[2] == "transposable_element" or lineAr[2] == "transposable_element_gene":
			if gffDict.get( chrm ) == None:
				gffDict[chrm] = {}
			gffDict[chrm][(start,end)] = strand
	
	gffFile.close()
	return gffDict

def processGFFDict( gffDict ):
	
	gffAr = []
	
	# loop through chromosomes
	for chrm in sorted(gffDict.keys()):
		# loop through positions - sorted
		print( '{:s}\t{:d}'.format( chrm, len( gffDict[chrm].keys() ) ) )
		for tup in sorted( gffDict[chrm].keys() ):
			start = tup[0]
			end = tup[1]
			strand = gffDict[chrm][tup]
			gffAr += [ (chrm, start, end, strand) ]
	return gffAr

def processFile( fileStr, outFileStr, gffAr, numBins, numBinsStream, streamSize, percentile, info ):
	print( 'Reading {:s}...'.format( fileStr ) )
	bedDict = readBed( fileStr )
	print( 'Processing {:s}'.format( fileStr ) )
	outMatrix, numList = processBed( bedDict, gffAr, numBins, numBinsStream, streamSize )
	thresh = calculatePercentile( percentile, numList )
	writeOutput( outFileStr, outMatrix, thresh, info )
	print( 'Output written for {:s}'.format( outFileStr ) )
	return 1

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
	
	bedFile.close()
	return bedDict

def processBed( bedDict, gffAr, numBins, numBinsStream, streamSize ):
	# fpkmDict organized {fpkm: (gene1, gene2, ...)}
	# geneDict organized {name: (scaffold, start, end, strand)}
	
	outMatrix = []
	numList = []
	# loop by gene length
	for info in gffAr:
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
		if curGeneVarAr == False:
			continue
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
		# add to numList
		numList += outAr
		outMatrix.append( outAr )
	return outMatrix, numList

def varByRegion( bedDict, chrm, start, end, numBins ):
	''' 
		takes in the variant dictionary generated for a gene and the start and 
		end positions of that gene. Note the dictionary is set up so the key is 	
		the position
		on the scaffold and the value is the frequency of that variant
		returns one array with number of variants separated by bin
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	if binWidth < 1:
		return False
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

	outFile = open( outFileStr, 'w' )
	geneCount = len( outMatrix )
	headerAr = ["te.num"] + ["T{:02}".format(i) for i in range(len(outMatrix[0])) ]
	outFile.write( info + "\n" + ",".join(headerAr) + "\n" )
	for gene in outMatrix:
		for i in range( 0, len(gene) ):
			if gene[i] > threshold:
				gene[i] = threshold
		outStr = 'TE{:08},'.format( geneCount )
		outStrAr = [ '{:.5f}'.format(x) for x in gene ]
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
	startInd = 0
	
	for i in range( min(4, len(argv) ) ):
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
	# end for
	
	gffFileStr = argv[startInd]
	bedFileStrAr = []
	
	for j in range( startInd + 1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 heatmap_te_bed_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
