import sys, math, glob, multiprocessing, subprocess, os, gzip

# Usage: python gbm_metaplot_pe.py [-q] [-h] [-m=meth_types] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-o=output_identifer] <gff_file> <allc_path> <sample1> [sampleN]*
## Calculate methylation level in specific sequence context for gene body methylation
## to produce metaplot(s)


NUMPROC=1
NUMBINS=20
STREAMSIZE=1000

def processInputs( gffFileStr, allcPath, sampleNamesAr, numProc, numBins, numBinsStream, streamSize, outID, methTypes, isPrint ):
	
	if isPrint:
		print( 'GFF File:', os.path.basename(gffFileStr) )
		print( 'AllC path:', allcPath )
		print( 'Samples:', ', '.join(sampleNamesAr) )
		print( 'Methylation types:',  ', '.join(methTypes) )
		print( 'Number of bins:', numBins )
		print( 'Number of bins stream:', numBinsStream )
		print( 'Stream size:', streamSize )
	
	if isPrint:
		print( 'Reading GFF ', os.path.basename(gffFileStr) )
		
	geneDict, cdsDict = readGFF( gffFileStr )

	if isPrint:
		print( 'Begin processing {:d} samples and {:d} chromosomes with {:d} processors...'.format( len(sampleNamesAr), len(list(geneDict.keys())), numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	
	results = [ pool.apply_async( processSample, args=(geneDict, cdsDict, allcPath, sample, numBins, numBinsStream, streamSize, methTypes, isPrint) ) for sample in sampleNamesAr ]
	outMatrix = [ p.get() for p in results ]
	info = '#from_script: gbm_metaplot_pe.py; numBins: {:d}; numBinsStream: {:d}; streamSize: {:d}; methTypes: {:s}'.format( numBins, numBinsStream, streamSize, ','.join(methTypes) )
	
	outFileStr = '{:s}_gbm_metaplot.tsv'.format(outID)
	
	if isPrint:
		print( 'Writing output to {:s}'.format( outFileStr) )
	writeOutput( outFileStr, outMatrix, info, sampleNamesAr, methTypes )
	print( 'Done.' )

def readGFF( gffFileStr ):
	'''
		Organized by chrm
		Return dictionary of gene information and CDS information
	'''
	geneDict = {}
	cdsDict = {}
	
	gffFile = open( gffFileStr, 'r' )
	
	# foundPrimary
	foundPrimary = False
	isPrimary = None
	
	for line in gffFile:
		lineAr = line.strip().split('\t')
		# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
		# (6) strand (7) frame (8) notes 
		if line.startswith('#') or len(lineAr) < 9:
			continue
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		strand = lineAr[6]
		feat = lineAr[2]
		notes = lineAr[8]
		
		if geneDict.get(chrm) == None:
			geneDict[chrm] = []
		if cdsDict.get(chrm) == None:
			cdsDict[chrm] = []
		
		# is gene
		if feat == 'gene':
			foundPrimary = False
			isPrimary = None
			geneDict[chrm] += [(start, end, strand)]
		# mrna
		elif feat == 'mRNA':
			# longest in notes -> one of them is labeled primary
			if 'longest=1' in notes:
				isPrimary = True
				foundPrimary = True
			elif 'longest=0' in notes:
				isPrimary = False
			# longest not labeled -> if isPrimary, unset and found primary	
			elif isPrimary:
				isPrimary = False
			# haven't found primary yet, this is primary
			elif foundPrimary == False:
				isPrimary = True
				foundPrimary = True
			# otherwise, we have primary already so keep going to the next gene
		# cds
		elif feat == 'CDS' and isPrimary:
			cdsDict[chrm] += [(start, end)]
	# end for line
	gffFile.close()
	return geneDict, cdsDict

def processSample( geneDict, cdsDict, allcPath, sampleName, numBins, numBinsStream, streamSize, methTypes, isPrint ):
	'''
		Produce values for a single sample
		returns Array of bins
	'''
	outAr = None
	chrmList = list(geneDict.keys())
	
	for chrm in chrmList:
		if isPrint:
			print('Processing {:s}: {:s}'.format( sampleName, chrm ) )
		allcFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format(allcPath, sampleName, chrm) )
		allcDict = readAllc( allcFileStr, cdsDict[chrm], methTypes )
		
		tmpAr = processChrm(geneDict[chrm], allcDict, numBins, numBinsStream, streamSize, methTypes)
		## add values to existing array
		outAr = addToArrayMeth( outAr, tmpAr, methTypes )
	# end for chrm
	
	#print( 'Processing {:s}...'.format( fileStr ) )
	#outMat = processAllC( gffDict, allCDict, numBins, numBinsStream, streamSize )
	return outAr

def readAllc( allcFileStr, cdsAr, methTypes ):
	
	outDict = {}
	
	if os.path.isfile( allcFileStr+'.gz' ):
		allcFile = gzip.open( allcFileStr + '.gz', 'rt')
	elif os.path.isfile( allcFileStr ):
		allcFile = open( allcFileStr, 'r' )
	else:
		print('ERROR:', allcFileStr, 'do not exist')
		return False
	
	cdsAr.sort()
	curIndex = 0
	curCds = cdsAr[curIndex]
	
	for line in allcFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		pos = int( lineAr[1] )
		
		# move to next cds
		while pos > curCds[1] and curIndex < (len(cdsAr)-1):
			curIndex += 1
			curCds = cdsAr[curIndex]
		
		# no more cds left in chrm
		if curIndex >= len(cdsAr):
			break
			
		# in the CDS
		if pos >= curCds[0] and pos <= curCds[1]:
			mType = findMethylType( lineAr[3] )
			mCount = int( lineAr[4] )
			tCount = int( lineAr[5] )
			if mType in methTypes or 'C' in methTypes:
				outDict[pos] = ( mCount, tCount, mType )
		
		# otherwise pos < curCds[0], move on
	# end for line
	allcFile.close()
	return outDict
	
def findMethylType( mc ):
	if mc.startswith( 'CG' ):
		return 'CG'
	elif mc.endswith( 'G' ):
		return 'CHG'
	elif mc == 'CNN':
		return 'CNN'
	else:
		return 'CHH'
		
def processChrm( geneAr, allcDict, numBins, numBinsStream, streamSize, methTypes ):
	
	totalGeneAr = emptyBinArray( numBins, methTypes )
	totalUpAr = []
	totalDownAr = []
	
	if STREAMSIZE != 0:
		totalUpAr = emptyBinArray( numBinsStream, methTypes )
		totalDownAr = emptyBinArray( numBinsStream, methTypes )
	
	# loop by gene
	for gene in geneAr:
		start = gene[0]
		end = gene[1]
		strand = gene[2]
		if STREAMSIZE != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = methByRegion(allcDict, upstream, start, numBinsStream, methTypes)
			# downstream
			curDownVarAr = methByRegion(allcDict, end, downstream, numBinsStream, methTypes)
		
		curGeneVarAr = methByRegion(allcDict, start, end, numBins, methTypes)
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr, methTypes )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr, methTypes )
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr, methTypes )
			else:
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr, methTypes )
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr, methTypes )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr, methTypes )
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr, methTypes )
			else:
				curGeneVarAr.reverse()
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr, methTypes )
	outAr = totalUpAr + totalGeneAr + totalDownAr
	return outAr

def methByRegion(allcDict, start, end, numBins, methTypes):
	''' 
		
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	outAr = []
	
	# loop for each bin
	for bin in range(numBins - 1):
		binDict = emptyBin( methTypes )
		# loop for each position in that bin
		for pos in range(binWidth):
			key = bin * binWidth + pos + start
			#print( key )
			try:
				dictEntry = allcDict.get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					mCount = dictEntry[0]
					tCount = dictEntry[1]
					mType = dictEntry[2]
					binDict[mType][0] += mCount
					binDict[mType][1] += tCount
					if 'C' in methTypes:
						binDict['C'][0] += mCount
						binDict['C'][1] += tCount
			except AttributeError:
				pass
		# end for pos
		outAr.append(binDict)
	# end for bin
	
	# handle the last "catch-all" bin
	lastBin = emptyBin( methTypes )
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = allcDict.get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				mCount = dictEntry[0]
				tCount = dictEntry[1]
				mType = dictEntry[2]
				lastBin[mType][0] += mCount
				lastBin[mType][1] += tCount
				if 'C' in methTypes:
					lastBin['C'][0] += mCount
					lastBin['C'][1] += tCount
		except AttributeError:
			pass
	# end for key
	outAr.append(lastBin)
	return outAr

def emptyBin( methTypes ):
	outDict = {}
	for m in methTypes:
		outDict[m] = [0, 0]
	return outDict

def emptyBinArray( nBins, methTypes ):
	return [emptyBin(methTypes) for x in range(nBins)]
			
def addToArrayMeth(oldAr, currAr, methTypes):
	if oldAr == None:
		return currAr
	if len(oldAr) != len(currAr):
		return -1
	for i in range( len(oldAr) ):
		#oldBinDict = oldAr[i]
		curBinDict = currAr[i]
		for m in methTypes:
			oldAr[i][m][0] += curBinDict[m][0]
			oldAr[i][m][1] += curBinDict[m][1]
		# end for m
	# end for i
	return oldAr

def writeOutput( outFileStr, outMatrix, info, sampleNamesAr, methTypes ):
	'''
		writes output to outFileStr
		output for all input files is written to one file
		array samples -> array bins -> dict methTypes -> tuple counts
	'''
	outFile = open( outFileStr, 'w' )
	
	# header
	outFile.write( info + '\nsample\tbin\tmethType\tmCount\ttCount\tvalue\n' )
	
	for j in range( len(sampleNamesAr) ):
		for i in range( len( outMatrix[j] ) ):
			binDict = outMatrix[j][i]
			for m in methTypes:
				mCount = binDict[m][0]
				tCount = binDict[m][1]
				wMeth = calcWeightedMeth( mCount, tCount )
				outStr = '{:s}\t{:d}\t{:s}\t{:d}\t{:d}\t{:f}\n'.format( sampleNamesAr[j], i, m, mCount, tCount, wMeth )
				outFile.write( outStr )
			# end for m
		# end for i
	# end for j
	outFile.close()

def calcWeightedMeth( mCount, tCount ):
	if mCount == 0 and tCount == 0:
		return float(0)
	elif tCount == 0:
		return float(-1) # error
	else:
		return float(mCount) / float(tCount)

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	methTypes = ['CG','CHG','CHH','C']
	isPrint = True
	outID = 'out'
	startInd = 0
	
	for i in range( min(6, len(argv) ) ):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-m='):
			methTypes = argv[i][3:].upper().split(',')
			for m in methTypes:
				if m not in ['CG', 'CHG', 'CHH', 'C']:
					print('ERROR: {:s} is not a valid methylation type'.format(m))
					exit()
			startInd += 1
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
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-h' :
			printHelp()
			exit()
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid argument'.format( argv[i] ) )
			exit()
	# end for
	
	gffFileStr = argv[startInd]
	allcPath = argv[startInd + 1]
	sampleNamesAr = []
	
	for j in range( startInd + 2, len(argv) ):
		sampleNamesAr += [ argv[j] ]
	
	processInputs( gffFileStr, allcPath, sampleNamesAr, numProc, numBins, numBinsStream, streamSize, outID, methTypes, isPrint )

def printHelp():
	print (" Usage: python3 gbm_metaplot_pe.py [-q] [-h] [-p=num_proc]  [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-o=output_identifer] [-m=meth_types] <gff_file> <allc_path> <sampleName1> [sampleNameN]*\n")
	print('Required:')
	print("gff_file\tpath to GFF formatted file\n\t\tmust have 'genes', 'mRNA', and 'CDS'" )
	print("bed_file\tpath to BED formatted file" )
	print( 'allc_file\tpath to allC file; all chrms in one file' )
	print( 'Optional:' )
	print( '-h\t\tprint this help menu and exit' )
	print( '-q\t\tquiet mode; do not print progress' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )
	print( '-b=num_bins[,num_bins_stream]\tnumber of bins to use for gene body and optionally up/down stream region [default 20]' )
	print( '-s=stream_size\tnumber of bp for up/down stream regions [default 1000]' )
	print( '-o=output_identifier\tstring included in output file name [default "out"]' )
	print("-m=meth_types\tcomma-separated list methylation contexts to analyze [default CG,CHG,CHH,C]")

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
