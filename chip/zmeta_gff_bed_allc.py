import sys, math, os, multiprocessing
NUMPROC=2

''' possible lines to change: 6, 9, 34, 188'''

# number of bins each section is broken into - change if desired
NUMBINS = 50
NUMBINS_STREAM=10

# size upstream and downstream of start/stop - change if desired
STREAMSIZE = 200

# usage: python3.4 meta_gff_bed.py <gff_file> <bed_file> [bed_file]*

def readGFF(gffFileStr):
	gffFile = open (gffFileStr, 'r' )
	geneAr = []
	
	for line in gffFile:
		line = line[:-1]
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		scaffold = lineAr[0]
		strand = lineAr[6]
		#print ( "GFF:",lineAr)
		
		# only apply to type = gene
		
		# change this to select a different gff feature 
		# or remove and un-indent line after
		if lineAr[2] == "gene":
			geneAr+= [(scaffold, start, end, strand)]
	
	gffFile.close()
	print( "Finished reading GFF file" )
	return geneAr
	
def processBed( geneAr, varDict, countReads ):
	# repeatAr organized [(scaffold, start, end, strand), ...]
	
	totalGeneAr = [0] * NUMBINS
	totalUpAr = []
	totalDownAr = []
	
	if STREAMSIZE != 0:
		totalUpAr = [0] * NUMBINS_STREAM
		totalDownAr = [0] * NUMBINS_STREAM
	
	# loop by gene length
	for gene in geneAr:
		#print( region )
		scaffold = gene[0]
		start = gene[1]
		end = gene[2]
		strand = gene[3]
		outAr = []
		if STREAMSIZE != 0:
			upstream = start - STREAMSIZE
			downstream = end + STREAMSIZE
			# upstream
			curUpVarAr = varByRegion(varDict, scaffold, upstream, start, NUMBINS_STREAM)
			# downstream
			curDownVarAr = varByRegion(varDict, scaffold, end, downstream, NUMBINS_STREAM)
		
		# repeat region
		curGeneVarAr = varByRegion(varDict, scaffold, start, end, NUMBINS)
		#print( 'curRepeatVarAr', curRepeatVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if STREAMSIZE != 0:
				#outAr = curUpVarAr + curRepeatVarAr + curDownVarAr
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				#outAr = curRepeatVarAr
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		else:
			if STREAMSIZE != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				#outAr = curDownVarAr + curRepeatVarAr + curUpVarAr
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalRepeatAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				curGeneVarAr.reverse()
				#outAr = curRepeatVarAr
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		#print( 'totalRepeatAr', totalRepeatAr )
	outAr = totalUpAr + totalGeneAr + totalDownAr
	outAr = [ float(x)/countReads*1000 for x in outAr ]
	return outAr

def readBed(bedFileStr ):
	''' takes in the bed file, scaffold we are considering
		returns a dictionary where the position is the key and the value is the number
		of hits for that position
	'''
	bedFile = open( bedFileStr, 'r' )
	varDict = {}
	countReads = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		line = line.replace('\n','')
		lineAr =line.split('\t')
		try:
			curScaf = lineAr[0]
			pos = int(lineAr[1])
		except ValueError:
			pass
		else:
			countReads += 1
			#print lineAr
			# no dictionary for scaffold
			if varDict.get(curScaf) == None:
				varDict[curScaf] = {pos:1}
			# dictionary for scaffold but position not included
			elif varDict.get(curScaf).get(pos) == None:
				varDict[curScaf][pos] = 1
			# dictionary for scaffold and position included
			else:
				varDict[curScaf][pos] += 1
			#print varDict[curScaf]
	
	bedFile.close()
	print( 'Finished reading', bedFileStr )
	return varDict, countReads

def varByRegion(varDict, scaffold, start, end, numBins):
	''' takes in the variant dictionary generated for a gene and the start and end
		positions of that gene. Note the dictionary is set up so the key is the position
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
				dictEntry = varDict.get(scaffold).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					varAr[bin] += dictEntry
			except AttributeError:
				#print "error with scaffold", scaffold
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = varDict.get(scaffold).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				varAr[-1] += dictEntry
		except AttributeError:
			#print "error with scaffold", scaffold
			pass
	return varAr

def methByRegion(allCDict, scaffold, start, end, numBins):
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
				dictEntry = allCDict.get(scaffold).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					methAr[bin]+= dictEntry[0]
					totalAr[bin] += dictEntry[1]
			except AttributeError:
				#print "error with scaffold", scaffold
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = allCDict.get(scaffold).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				methAr[-1] += dictEntry[0]
				totalAr[-1] += dictEntry[1]
		except AttributeError:
			#print ( "error with scaffold", scaffold )
			pass
	z = zip( methAr, totalAr )
	outAr = [ (m,t) for m,t in z ]
	#print (outAr)
	return outAr

def readAllC( allCFileStr ):
	'''
		This allC file should be the allC information for all chromosomes in one file
		not just a single chromosome
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
			if chrm.startswith('scaffold')==False:
				chrm = 'Chr{:02d}'.format( int( lineAr[0] ) )
			if allCDict.get( chrm ) == None:
				allCDict[ chrm ] = {}
			allCDict[chrm][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
	
	allCFile.close()
	print( len( allCDict['Chr01'] ) )
	print( 'Finished reading', allCFileStr )
	return allCDict

def processAllC( geneAr, allCDict ):
	
	totalGeneAr = [(0,0)] * NUMBINS
	totalUpAr = []
	totalDownAr = []
	
	if STREAMSIZE != 0:
		totalUpAr = [(0,0)] * NUMBINS_STREAM
		totalDownAr = [(0,0)] * NUMBINS_STREAM
	
	# loop by gene length
	for gene in geneAr:
		#print( region )
		scaffold = gene[0]
		start = gene[1]
		end = gene[2]
		strand = gene[3]
		outAr = []
		if STREAMSIZE != 0:
			upstream = start - STREAMSIZE
			downstream = end + STREAMSIZE
			# upstream
			curUpVarAr = methByRegion(allCDict, scaffold, upstream, start, NUMBINS_STREAM)
			# downstream
			curDownVarAr = methByRegion(allCDict, scaffold, end, downstream, NUMBINS_STREAM)
		
		curGeneVarAr = methByRegion(allCDict, scaffold, start, end, NUMBINS)
		#print( 'curRepeatVarAr', curRepeatVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if STREAMSIZE != 0:
				#outAr = curUpVarAr + curRepeatVarAr + curDownVarAr
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr )
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
			else:
				#outAr = curRepeatVarAr
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
		else:
			if STREAMSIZE != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				#outAr = curDownVarAr + curRepeatVarAr + curUpVarAr
				totalUpAr = addToArrayMeth( totalUpAr, curUpVarAr )
				totalDownAr = addToArrayMeth( totalDownAr, curDownVarAr )
				totalRepeatAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
			else:
				curGeneVarAr.reverse()
				#outAr = curRepeatVarAr
				totalGeneAr = addToArrayMeth( totalGeneAr, curGeneVarAr )
	outAr = totalUpAr + totalGeneAr + totalDownAr
	#print( outAr )
	methAr = calculateMethylation( outAr )
	return methAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	
	return oldAr

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
	
	outAr = []
	for i in range(len(methAr)) :
		outAr += [ (float(methAr[i][0])/float(methAr[i][1]))*100]
	return outAr

def writeOutput( varMatrix, bedSampleNames ):
	'''
		writes output to predetermined file name
		output for all input bed files is written to one file
	'''
	# can change outfile name or outfile prefix
	outFileStr = "meta_"+str(NUMBINS)+"_"+str(STREAMSIZE)+".csv"
	
	outFile = open( outFileStr, 'w' )
	
	# header
	outFile.write( 'bin,count,type\n' )
	
	# for each bed file
	for j in range( len(bedSampleNames) ):
		
		for i in range( len( varMatrix[j] ) ):
			outStr = "{:d},{:.5f},{:s}\n".format( i+1, varMatrix[j][i], bedSampleNames[j] )
			outFile.write( outStr )
	print( "Output written" )
	outFile.close()

def checkFile( fileStr ):
	
	if fileStr.endswith( '.bed' ):
		return True
	else:
		# allC file
		return False

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

def processFile( gffAr, fileStr ):
	isBed = checkFile( fileStr )
	if isBed:
		bedDict, countReads = readBed( fileStr )
		o = processBed( gffAr, bedDict, countReads )
	else:
		allCDict = readAllC( fileStr )
		o = processAllC( gffAr, allCDict )
	print( 'Finished processing', fileStr )
	return o

def processInput( gffFileStr, fileStrAr ):
	''' 
		calls the readRepeat to get the define regions of repeats
		then processes those regions onto the bed files
		then calls to write with appropriate header
	'''
	print( str(NUMBINS), 'number of bins per section' )
	print( str(STREAMSIZE), 'upsteam and downstream size' )
	sampleNames = getSampleNames( fileStrAr )
	
	gffAr = readGFF( gffFileStr )
	
	pool = multiprocessing.Pool( processes=NUMPROC )
	print( 'Begin processing files with {:d} processors'.format( NUMPROC ) )
	results = [ pool.apply_async( processFile, args=(gffAr, f) ) for f in fileStrAr ]
	varMatrix = [ p.get() for p in results ]
	
	writeOutput( varMatrix, sampleNames )
		#cleanOutputFile( outFileStr )
	
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print( "Usage: meta_gff_bed_allC.py [-m=meth_types] <gff_file> <bed_file> [bed_file]*" )
	else:
		gffFileStr = sys.argv[1]
		n = len(sys.argv) - 2
		bedFileStrAr = [None] * n
		# get the bed file names
		for s in range (2, len(sys.argv)):
			#print (s,sys.argv[s])
			bedFileStrAr[s-2] = sys.argv[s]
		processInput(gffFileStr, bedFileStrAr)
