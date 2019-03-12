import sys, math, glob, multiprocessing, subprocess, os
from bioFiles import *

# Usage: python3 heatmap_bed_allc_fr_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*

NUMPROC=1
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95

def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile ):

	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nPercentile for correction: {:.3f}\nSamples included: {:s}\n'.format( numBins, numBinsStream, streamSize, percentile, ', '.join(sampleNamesAr) ) )
	print( 'Reading FPKM file' )
	fpkmFile = FileFPKM( fpkmFileStr )
	fpkmAr = fpkmFile.getFPKMArray( )
	
	print( 'Reading GFF file' )
	gffFile = FileGFF( gffFileStr )
	geneDict, chrmFormat  = gffFile.getGeneDict( chrm=True )
	print( chrmFormat )
	outFileStrAr = [ '{:s}_{:d}_{:d}_{:d}.csv'.format( sampleName, streamSize, numBins, int(percentile*100) ) for sampleName in sampleNamesAr ]
	
	info = "#from_script:heatmap_bed_allc_fr_pe.py; num_bins:{:d}; num_bins_stream:{:d}; stream_size:{:d}; percentile:{:.2f}; gff_file:{:s}; fpkm_file:{:s}".format( numBins, numBinsStream, streamSize, percentile, str(gffFile), str(fpkmFile) )
	
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	results = [ pool.apply_async( processFile, args=(bedFileStrAr[i], outFileStrAr[i], fpkmAr, geneDict, numBins, numBinsStream, streamSize, percentile, chrmFormat, info )) for i in range( len(bedFileStrAr) ) ]
	fin = [ p.get() for p in results ]
	if sum(fin) != len( bedFileStrAr ):
		print( 'ERROR: not enough files written' )
		exit()
	print( 'Done' )

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		if fileStr.endswith( '.bed' ):
			leftIndex = fileStr.rfind('/')
			rightIndex = fileStr.rfind('.')
			sampleName = fileStr[leftIndex+1:rightIndex]
			sampleNamesAr += [ sampleName ]
		else:
			# allC file
			leftIndex = fileStr.rfind('/')
			sampleName = fileStr[leftIndex+1:]
			aIndex = sampleName.find( '_' )
			bIndex = sampleName[aIndex+1:].find( '_' )
			sampleNamesAr += [ sampleName[:bIndex+aIndex+1]+'_mCG' ]
	return sampleNamesAr

def processFile( fileStr, outFileStr, fpkmAr, gffDict, numBins, numBinsStream, streamSize, percentile, chrmFormat, info ):
	isBed = checkFile( fileStr )
	if isBed:
		print( 'Processing {:s}'.format( fileStr ) )
		bedFile = FileBED_FR( fileStr )
		bedDict, bpCount = bedFile.getBedDict( )
		outMatrix, numList = processBed( bedDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize, bpCount )
		thresh = calculatePercentile( percentile, numList )
		writeOutput( outFileStr, outMatrix, thresh, info )
		print( 'Output written for {:s}'.format( outFileStr ) )
		return 1
	else:
		print( 'Reading {:s}...'.format( fileStr ) )
		allCDict = readAllC( fileStr, chrmFormat )
		for key in allCDict.keys():
			print( key, len(allCDict[key]) )
		print( 'Processing {:s}'.format( fileStr ) )
		outMatrix, numList = processAllC( allCDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize )
		thresh = calculatePercentile( percentile, numList )
		writeOutput( outFileStr, outMatrix, thresh, info )
		print( 'Output written for {:s}'.format( outFileStr ) )
		return 1

def checkFile( fileStr ):
	if fileStr.endswith( '.bed' ) or fileStr.endswith( '.bed.gz' ):
		return True
	else:
		# allC file
		return False

def processBed( bedDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize, bpCount ):
	# fpkmDict organized {fpkm: (gene1, gene2, ...)}
	# geneDict organized {name: (scaffold, start, end, strand)}
	
	outMatrix = []
	numList = []
	# loop by gene length
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
		# adjust by bpCount
		outAr = [ x / bpCount * 1000000 for x in outAr ]
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

def readAllC( allCFileStr, chrmFormat ):
	'''
		This allC file should be the allC information for all chromosomes in one 
		file not just a single chromosome
	'''
	allCFile = open( allCFileStr, 'r' )
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

def processAllC( allCDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize ):
	
	outMatrix = []
	numList = []
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
		numList += outAr
	return outMatrix, numList

def methByRegion(allCDict, chrm, start, end, numBins):
	''' 
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print ('binWidth', binWidth)
	methAr = [0] * numBins
	totalAr = [0] * numBins
	
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

def calculateMethylation( methAr, length ):
	
	outAr = [-1] * len( methAr )
	for i in range(len(methAr)) :
		if methAr[i][0] == 0 and methAr[i][1] == 0:
			outAr[i] = 0
		elif methAr[i][0] != 0 and methAr[i][1] == 0:
			print( methAr[i] )
		else:
			outAr[i] = (float(methAr[i][0]) / float(methAr[i][1]))* 100 / length
	return outAr

def writeOutput( outFileStr, outMatrix, threshold, info ):

	outFile = open( outFileStr, 'w' )
	geneCount = len( outMatrix )
	headerAr = ["gene.num"] + ["B{:02}".format(i) for i in range(len(outMatrix[0])) ]
	outFile.write( info + "\n" + ",".join(headerAr) + "\n" )
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
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3 heatmap_bed_allc_fr_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*")
	else:
		parseInputs( sys.argv[1:] )
