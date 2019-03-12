import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 metaplot_from_bed_allc_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-o=output_identifer] [-te] <gff_file> <bed_file> [bed_file | allc_file]*


NUMPROC=2
NUMBINS=20
STREAMSIZE=1000

def processInputs( gffFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, outID, isTE ):

	sampleNamesAr = getSampleNames( bedFileStrAr )
	if outID == '' and isTE == False:
		outFileStr = 'meta_{:d}_{:d}.csv'.format(numBins,streamSize)
	elif outID == '' and isTE:
		outFileStr = 'meta_te_{:d}_{:d}.csv'.format(numBins,streamSize)
	elif isTE:
		outFileStr = 'meta_te_{:d}_{:d}_{:s}.csv'.format( numBins,streamSize,outID )
	else:
		outFileStr = 'meta_{:d}_{:d}_{:s}.csv'.format(numBins,streamSize,outID)
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nFeatures: {:s}\nSamples included: {:s}\nOutput written to: {:s}'.format( numBins, numBinsStream, streamSize,('transposons' if isTE else 'genes'), ', '.join(sampleNamesAr), outFileStr ) )
	
	print( 'Reading GFF file {:s}...'.format( gffFileStr ) )
	gffAr, chrmFormat = readGFF( gffFileStr, isTE )
	
	print( 'Begin processing files with {:d} processors...'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(gffAr, f, numBins, numBinsStream, streamSize, chrmFormat) ) for f in bedFileStrAr ]
	varMatrix = [ p.get() for p in results ]
	
	print( 'Writing output...' )
	writeOutput( outFileStr, varMatrix, sampleNamesAr )
	print( 'Done.' )

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
			sampleNamesAr += [ 'mCG methylation' ]
	return sampleNamesAr

def readGFF(gffFileStr, isTE):
	gffFile = open (gffFileStr, 'r' )
	geneAr = []
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
		#print ( "GFF:",lineAr)
		
		# only apply to type = gene
		
		# change this to select a different gff feature 
		# or remove and un-indent line after
		if isTE and lineAr[2] in ['transposable_element', 'transposable_element_gene']:
			geneAr+= [(chrm, start, end, strand)]
		elif isTE == False and lineAr[2] == "gene":
			geneAr+= [(chrm, start, end, strand)]
	
	gffFile.close()
	print( "Finished reading GFF file" )
	return geneAr, chrmFormat

def processFile( gffAr, fileStr, numBins, numBinsStream, streamSize, chrmFormat ):
	isBed = checkFile( fileStr )
	if isBed:
		print( 'Reading BED file {:s}...'.format( fileStr ) )
		bedDict, countReads = readBed( fileStr )
		print( 'Processing {:s}...'.format( fileStr ) )
		outMat = processBed( gffAr, bedDict, countReads, numBins, numBinsStream, streamSize )
	else:
		print( 'Reading allC file {:s}...'.format( fileStr ) )
		allCDict = readAllC( fileStr, chrmFormat )
		print( 'Processing {:s}...'.format( fileStr ) )
		outMat = processAllC( gffAr, allCDict, numBins, numBinsStream, streamSize )
	print( 'Finished processing', fileStr )
	return outMat

def checkFile( fileStr ):
	
	if fileStr.endswith( '.bed' ):
		return True
	else:
		# allC file
		return False

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

def processBed( geneAr, bedDict, countReads, numBins, numBinsStream, streamSize ):
	# repeatAr organized [(scaffold, start, end, strand), ...]
	
	totalGeneAr = [0] * numBins
	totalUpAr = []
	totalDownAr = []
	
	if streamSize != 0:
		totalUpAr = [0] * numBinsStream
		totalDownAr = [0] * numBinsStream
	
	# loop by gene length
	for gene in geneAr:
		#print( region )
		chrm = gene[0]
		start = gene[1]
		end = gene[2]
		strand = gene[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = varByRegion(bedDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = varByRegion(bedDict, chrm, end, downstream, numBinsStream)
		
		# gene region
		curGeneVarAr = varByRegion(bedDict, chrm, start, end, numBins)
		if curGeneVarAr == False:
			continue
		#print( curGeneVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				#outAr = curDownVarAr + curRepeatVarAr + curUpVarAr
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				curGeneVarAr.reverse()
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		# end if-else strand
	# end for gene
	outAr = totalUpAr + totalGeneAr + totalDownAr
	#print( outAr )
	outAr = [ float(x)/countReads*1000 for x in outAr ]
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
	return varAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	return oldAr

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
		if lineAr[3].startswith( 'CG' ) and lineAr[6] == '1':
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

def processAllC( geneAr, allCDict, numBins, numBinsStream, streamSize ):
	
	totalGeneAr = [(0,0)] * numBins
	totalUpAr = []
	totalDownAr = []
	
	if STREAMSIZE != 0:
		totalUpAr = [(0,0)] * numBinsStream
		totalDownAr = [(0,0)] * numBinsStream
	
	# loop by gene length
	for gene in geneAr:
		#print( region )
		chrm = gene[0]
		start = gene[1]
		end = gene[2]
		strand = gene[3]
		outAr = []
		if STREAMSIZE != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = methByRegion(allCDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = methByRegion(allCDict, chrm, end, downstream, numBinsStream)
		
		curGeneVarAr = methByRegion(allCDict, chrm, start, end, numBins)
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr )
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
			else:
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr )
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
			else:
				curGeneVarAr.reverse()
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
	outAr = totalUpAr + totalGeneAr + totalDownAr
	methAr = calculateMethylation( outAr )
	return methAr

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
	outAr = [ (m,t) for m,t in z ]
	return outAr

def addToArrayMeth(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			mCount = oldAr[i][0] + currAr[i][0]
			tCount = oldAr[i][1] + currAr[i][1]
			oldAr[i] = (mCount, tCount)
	return oldAr
	
def calculateMethylation( methAr ):
	
	outAr = [-1] * len( methAr )
	for i in range(len(methAr)) :
		if methAr[i][0] == 0 and methAr[i][1] == 0:
			outAr[i] = 0
		elif methAr[i][0] != 0 and methAr[i][1] == 0:
			print( methAr[i] )
		else:
			outAr[i] = float(methAr[i][0]) / float(methAr[i][1]) * 100
	return outAr

def writeOutput( outFileStr, varMatrix, sampleNamesAr ):
	'''
		writes output to outFileStr
		output for all input files is written to one file
	'''
	outFile = open( outFileStr, 'w' )
	
	# header
	outFile.write( 'bin,count,sample\n' )
	
	# for each bed file
	for j in range( len(sampleNamesAr) ):
		
		for i in range( len( varMatrix[j] ) ):
			outStr = "{:d},{:f},{:s}\n".format( i+1, varMatrix[j][i], sampleNamesAr[j] )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	isTE = False
	outID = ''
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
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-te':
			isTE = True
			startInd += 1
	# end for
	
	gffFileStr = argv[startInd]
	bedFileStrAr = []
	
	for j in range( startInd + 1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, outID, isTE )
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print (" Usage: python3.4 metaplot_from_bed_allc_pe.py [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-o=output_identifer] [-te] <gff_file> <bed_file> [bed_file | allc_file]*")
	else:
		parseInputs( sys.argv[1:] )
