import sys, math, glob, multiprocessing, subprocess, os, bisect, random, pickle

# Usage: python3 metaplot_abs_chip_fr_pe.py [-p=num_proc] [-b=num_bins | -s=bin_size] [-d=distance] [-o=out_id] <gff_file> <bed_file> [bed_file]*
# using full read to generate absolute distance based metaplot
NUMPROC=1
NUMBINS=20
DISTANCE=1000
OUTID=''

def processInputs( gffFileStr, bedFileStrAr, numBins, binSize, numProc, outID, distance ):
	sampleNamesAr = getSampleNames( bedFileStrAr )
	
	#info
	info = '#from_script: metaplot_abs_chip_fr_pe.py; '
	if numBins == -1:
		info += 'bin_size: {:d}; '.format( binSize )
		numBins = int( float(distance) / binSize )
	else:
		info += 'num_bins: {:d}; '.format( numBins )
		binSize = int( float(distance) / numBins )
	info += 'distance: {:d}; gff_file: {:s}; samples: {:s}'.format( distance, os.path.basename(gffFileStr), ','.join(sampleNamesAr) )
	
	if outID == '':
		outFileStr = 'meta_abs_n{:d}_d{:d}.tsv'.format( numBins, distance )
	else:
		outFileStr = 'meta_abs_{:s}.tsv'.format( outID )
		
	print( 'GFF file:\t{:s}\nDistance:\t{:d}\nNum Bins:\t{:d}\nBin size:\t{:d}\nSamples:\t{:s}\nOut file:\t{:s}\n'.format( os.path.basename(gffFileStr), distance, numBins, binSize, outFileStr, ', '.join(sampleNamesAr) ))
	
	print( 'Reading GFF' )
	gffAr = readGFF( gffFileStr )
	
	# parallel processing
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(bedFile, gffAr, distance, numBins, binSize) ) for bedFile in bedFileStrAr ]
	outMatrix = [ p.get() for p in results ]
	
	# write ouput
	print( 'Writing output to {:s}'.format( outFileStr ) )
	writeOutput( outFileStr, outMatrix, info, sampleNamesAr )
	print( 'Done' )
	
def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr
	
def readGFF( inFileStr ):
	gffFile = open (inFileStr, 'r' )
	gffAr = []
	
	for line in gffFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		chrm = lineAr[0]
		strand = lineAr[6]
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			gffAr += [(chrm, start, end, strand)]
	return gffAr

def processFile( inFileStr, gffAr, distance, numBins, binSize ):
	print( 'Processing {:s}'.format(os.path.basename(inFileStr)) )
	# read BED file
	bedDict, bpCount = checkPickleFile( inFileStr )
	print( 'Analyzing {:s}'.format(os.path.basename(inFileStr)) )
	valueAr = processBed( bedDict, gffAr, distance, numBins, binSize, bpCount )
	return valueAr

def checkPickleFile( inFileStr ):
	'''
		check file for pickled version
		if pickled version exists, check modification time compared to file
		if picked version doesn't exist or out of date, create picked version
		note: to properly store dictionary with bp counts, the key '##counts##'
		will be added to the dictionary
		return bedDict and bpCount
	'''
	rInd = inFileStr.rfind( '.' )
	pickleFileStr = inFileStr[:rInd] + '.pick'
	
	# check existence
	fileExists = os.path.isfile( pickleFileStr )
	# if exists -> check modification time
	if fileExists:
		mTime = os.path.getmtime( inFileStr )
		tTime = os.path.getmtime( pickleFileStr )
		fSize = os.path.getsize( pickleFileStr )
		if tTime < mTime or fSize < 100:
			print( 'Updating {:s}...'.format(os.path.basename(pickleFileStr)) )
			bedDict = createPickleFile( inFileStr, pickleFileStr )
		else:
			print( 'Loading {:s}...'.format(os.path.basename(pickleFileStr)) )
			bedDict = pickle.load( open(pickleFileStr, 'rb') )
	# does not exist -> create
	else:
		print( 'Creating {:s}...'.format(os.path.basename(pickleFileStr)) )
		bedDict = createPickleFile( inFileStr, pickleFileStr )
	# get bp count
	counts = bedDict.get( '##counts##' )
	return bedDict, counts

def createPickleFile( inFileStr, pickleFileStr ):
	# create dictionary
	print( 'Reading {:s}...'.format( os.path.basename(inFileStr) ) )
	bedDict, bpCounts = readBedFR( inFileStr )
	# add count to dictionary as '##counts##'
	bedDict[ '##counts##' ] = bpCounts
	# pickle
	pickle.dump( bedDict, open( pickleFileStr, 'wb' ), protocol=4 )
	return bedDict

def readBedFR( bedFileStr ):
	'''
		creates a dictionary with each scaffold and those dictionaries are a 
		dictionary
		for the positions of coordinates with frequency count from the bed file
		{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
		return the dictionary
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	count = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		lineAr =line.rstrip().split('\t')
		for pos in range( int(lineAr[1])+1, int(lineAr[2])+2):
			try:
				curScaf = lineAr[0]
				pos = int( lineAr[1] ) + 1
			except ValueError:
				pass
			else:
				count += 1
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
	return bedDict, count

def processBed( bedDict, gffAr, distance, numBins, binSize, bpCount ):
	tssAr = [0] * (2*numBins)
	ttsAr = [0] * (2*numBins)
	
	# loop through genes
	for i in range(len(gffAr)):
		chrm, start, end, strand = gffAr[i]
		chrmDict = bedDict.get(chrm)
		if (end-start+1) < distance or chrmDict == None:
			pass
		
		# tss
		tssStart = start - distance
		tssEnd = start + distance
		tss = countRegion( chrmDict, tssStart, start, tssEnd, numBins, binSize )
		
		# tts
		ttsStart = end - distance
		ttsEnd = end + distance
		tts = countRegion( chrmDict, ttsStart, end, ttsEnd, numBins, binSize )
		#if tss == -1 or tts == -1:
		#print( 'error with', chrm, start, end )
		try:
			u = tss[1]
			u = tts[2]
		except TypeError:
			print( 'error with int:', tss, tts, 'chrm', chrm )
		else:
			# handle negative strand genes
			if strand == '-':
				tss.reverse()
				tts.reverse()
				tssAr = addToArray( tssAr, tts )
				ttsAr = addToArray( ttsAr, tss )
			# positive genes
			else:
				tssAr = addToArray( tssAr, tss )
				ttsAr = addToArray( ttsAr, tts )
	# adjust by number of bp in library	(and bin size?)
	outAr = tssAr + ttsAr
	outs = [ float(x) / bpCount * 1000 for x in outAr ]
	return outs

def countRegion( bedDict, start, mid, stop, numBins, binWidth ):
	
	upAr = [0] * numBins
	downAr = [0] * numBins
	if bedDict == None:
		return -1
	# upstream Ar
	for bin in range( numBins ):
		bStart = int(start + bin*binWidth)
		# loop through each position
		for pos in range(bStart, int(bStart+binWidth)):
			dictEntry = bedDict.get( pos )
			if dictEntry != None:
				upAr[bin] += dictEntry
	
	# downstream Ar
	for bin in range(numBins):
		bStart = int(mid+1 + bin*binWidth)
		for pos in range(bStart, int(bStart+binWidth)):
			dictEntry = bedDict.get( pos )
			if dictEntry != None:
				downAr[bin] += dictEntry
	outAr = upAr + downAr
	return outAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	return oldAr

def writeOutput( outFileStr, outMatrix, info, sampleNamesAr ):
	
	outFile = open( outFileStr, 'w' )
	header = '\nsample\tregion\tbin\tvalue\n'
	outFile.write( info + header )
	
	# loop through matrix -> samples
	for i in range(len(outMatrix)):
		sampleAr = outMatrix[i]
		nBins = int(len(sampleAr) / 2)
		if nBins % 2 != 0:
			print ( 'error: nbins/2' )
		# loop through tss bins
		for j in range(nBins):
			outStr = '{:s}\t{:s}\t{:d}\t{:.5f}\n'.format( sampleNamesAr[i], 'tss', j, sampleAr[j] )
			outFile.write( outStr )
		# add intermediate rows
		outStr = '{0:s}\ttss\t{1:d}\tNA\n{0:s}\ttts\t{2:d}\tNA\n'.format( sampleNamesAr[i], nBins, nBins+1 )
		outFile.write( outStr )
		# loop through tts bins
		for j in range(nBins):
			outStr = '{:s}\t{:s}\t{:d}\t{:.5f}\n'.format( sampleNamesAr[i], 'tts', j+nBins+2, sampleAr[j+nBins] )
			outFile.write( outStr )
	# end loop

def parseInputs( argv ):
	# Usage: python3 metaplot_abs_chip_fr_pe.py [-p=num_proc] [-b=num_bins | -s=bin_size] [-d=distance] [-o=out_id] <gff_file> <bed_file> [bed_file]*
	
	numProc = NUMPROC
	binSize = -1
	numBins = -1
	distance = DISTANCE
	outID = OUTID
	startInd = 0
	
	for i in range(min(4,len(argv)-2)):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i].startswith( '-b=' ):
			# check for previous '-s'
			if binSize != -1:
				print( 'ERROR: cannot specify -b and -s together' )
				exit()
			try:
				numBins = int( argv[i][3:] )
				startInd += 1
				isNumBins = True
			except ValueError:
				print( 'ERROR: number of bins must be an integer' )
				exit()
		elif argv[i].startswith( '-s=' ):
			# check previous '-b'
			if numBins != -1:
				print( 'ERROR: cannot specify -b and -s together' )
				exit()
			try:
				binStr = argv[i][3:]
				if binStr.endswith( 'k' ) or binStr.endswith( 'K' ):
					binSize = int( binStr[:-1] ) * 1000
				elif binStr.endswith( 'm' ) or binStr.endswith( 'M' ):
					binSize = int( binStr[:-1] ) * 1000000
				else:
					binSize = int( binStr )
				startInd += 1
			except ValueError:
				print( 'ERROR: bin size must be an integer' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-d=' ):
			try:
				distance = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: distance must be an integer' )
				exit()
		elif argv[i].startswith('-'):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	
	if numBins == -1 and distance % binSize != 0:
		nnbins = int( math.ceil( float(distance) / binSize ) )
		distance = nnbins * binSize
		print( 'WARNING: bin size doesn\'t evenly divide into distance...adjusting distance to {:s}'.format(distance) )
	elif binSize == -1 and distance % numBins != 0:
		nnsize = int( math.ceil( float(distance) / numBins ) )
		distance = nnsize * numBins
		print( 'WARNING: num bins doesn\'t evenly divide into distance...adjusting distance to {:s}'.format(distance) )
	
	gffFileStr = argv[startInd]
	bedFileStrAr = []
	
	for j in range( startInd + 1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	processInputs( gffFileStr, bedFileStrAr, numBins, binSize, numProc, outID, distance )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3 metaplot_abs_chip_fr_pe.py [-p=num_proc] [-b=num_bins | -s=bin_size] [-d=distance] [-o=out_id] <gff_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
