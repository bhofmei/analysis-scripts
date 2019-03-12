import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: chrom_plot_bed_combine_pe.py [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-t=percentile] <fasta_index> <bed_file> [bed_file]*

NUMPROC=2
NUMBINS=100
PERCENTILE=0.95

def processInputs( fastaIndexStr, bedFileStrAr, numBins, binSize, chrmList, numProc, outID, percentile ):
	print( 'Reading FASTA index.' )
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	#print( chrmDict )
	print( 'Number of chromosomes: ' + str(len(chrmDict.keys())) )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	#print( aNumBins )
	#print( chrmDict )
	info = "#from_script:chrom_plot_bed_combine_pe.py;"
	if(binSize == -1):
		info += "num_bins:{:d}".format( aNumBins )
	else:
		info += "bin_size:{:d}".format( binSize )
	info += "; num_chrms:{:d}; percentile:{:.1f}".format( len( chrmDict.keys() ), percentile * 100 )
	outFileStr = 'chrm_bed'
	if outID != '':
		outFileStr += '_' + outID
	if numBins != -1:
		outFileStr += '_{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:d}'.format( binSize )
	outFileStr += '_{:d}'.format( int( percentile * 100 ) )
	outFileStr += '.tsv'
	
	# get individual matrices
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( readBedFile, args=(f, chrmDict, aNumBins, binSize ) ) for f in bedFileStrAr ]
	bedDictAr = [ p.get() for p in results ]
	
	# combine matrices
	bedDict = combineBedDict( bedDictAr )
	
	# determine percentile
	thresh = determinePercentile( bedDict, percentile )
	
	# output
	bedStr = bedDictToStr( bedDict, thresh )

	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, bedStr, info )
	
	print( 'Done.' )

def readFastaIndex( fastaIndexStr, chrmList ):
	
	indexFile = open( fastaIndexStr, 'r' )
	chrmDict = {}
	
	for line in indexFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) length ...
		chrm = lineAr[0]
		length = int( lineAr[1] )
		if chrmList == None or chrm in chrmList:
			chrmDict[chrm] = length
	indexFile.close()
	return chrmDict

def determineNumBins( chrmDict, numBins, binSize ):
	
	if binSize == -1:
		# correct for bin size
		for chrm in chrmDict.keys():
			aBinSize = math.ceil(chrmDict[chrm] / numBins)
			chrmDict[chrm] = aBinSize
		return numBins, chrmDict
	
	maxLength = 0
	for chrm in chrmDict.keys():
		if chrmDict[chrm] > maxLength:
			maxLength = chrmDict[chrm]
	
	n = int( math.ceil( float(maxLength) / binSize ) )
	return n, chrmDict

def readBedFile( bedFileStr, chrmDict, numBins, binSize ):
	
	bedDict = {}
	countReads = 0
	
	# add chrms to bedDict - each chrm is array (length based on numBins/binSize)
	for chrm in chrmDict.keys():
		if binSize == -1:
			bedDict[chrm] = [0] * numBins
		else:
			n = math.ceil( chrmDict[chrm] / binSize )
			bedDict[chrm] = [0]*n + [-1]*(numBins - n)
	
	print( 'Reading {:s}...'.format( bedFileStr ) )
	bedFile = open( bedFileStr, 'r' )
	
	for line in bedFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
		chrm = lineAr[0]
		pos = int( lineAr[1] ) + 1
		# want to skip chromosomes we have no info for
		if bedDict.get( chrm ) != None:
			if binSize == -1:
				bin = int( pos // chrmDict[chrm] )
			else:
				bin = int( pos // binSize )
			bedDict[chrm][bin] += 1
			countReads += 1
			
	bedFile.close()
	# adjust by number of reads
	#bedDict = adjustBedDict( bedDict, countReads )
	# adjust by library size and bin width
	bedDict = adjustBed( bedDict, countReads, binSize, chrmDict )
	return bedDict

def adjustBedDict( bedDict, readCounts ):
	for chrm in bedDict.keys():
		getAr = bedDict[ chrm ]
		newAr = [ ( -1 if x == -1 else float(x) / readCounts * 1000000 ) for x in getAr ]
		bedDict[chrm] = newAr
	return bedDict

def adjustBed( bedDict, readCounts, binSize, chrmDict ):
	for chrm in bedDict.keys():
		getAr = bedDict[ chrm ]
		if binSize == -1:
			divSize = chrmDict[ chrm ]
		else:
			divSize = binSize
		newAr = [ ( -1 if x == -1 else float(x) /( readCounts * divSize) * 10**9 ) for x in getAr ]
		bedDict[ chrm ] = newAr
	return bedDict
def combineBedDict( bedDictAr ):
	outDict = {}
	n = len( bedDictAr )
	#print( n)
	for chrm in bedDictAr[0].keys():
		#print( chrm, len( bedDictAr[0][chrm] ))
		newAr = bedDictAr[0][chrm]
		# sum
		for k in range( 1, n ):
			kAr = bedDictAr[k][chrm]
			newAr = [ ( -1 if newAr[i] == -1 else newAr[i] +kAr[i] ) for i in range(len(kAr)) ]
		#print( chrm, newAr[-10:] )
		# divide
		newAr2 = [ (-1 if x ==-1 else (float(x) / n) ) for x in newAr ]
		outDict[ chrm ] = newAr2
	# end for chrm
	return outDict

def determinePercentile( bedDict, percentile ):
	
	numList = []
	
	for chrm in bedDict.keys():
		tmpAr = bedDict[chrm]
		for i in range(len(tmpAr)):
			if tmpAr[i] != -1:
				numList += [tmpAr[i]]
	# end for
	
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
	
def bedDictToStr( bedDict, thresh ):
	
	outStr = ''
	# loop through chroms
	for chrm in sorted( bedDict.keys() ):
		tmpAr = bedDict[chrm]
		# loop through bins
		for i in range( len(tmpAr) ):
			if tmpAr[i] == -1:
				outStr += '{:s}\t{:d}\t\tNA\n'.format( chrm, i )
			elif tmpAr[i] > thresh:
				outStr += '{:s}\t{:d}\t{:.4f}\n'.format( chrm, i, float(thresh) )
			else:
				outStr += '{:s}\t{:d}\t{:.4f}\n'.format( chrm, i, float(tmpAr[i]) )
	return outStr

def writeOutput( outFileStr, outStr, info ):
	
	# header would be: chrm, bin, value
	
	outFile = open( outFileStr, 'w' )
	# header
	header = "#chrm\tbin\tvalue\n"
	outFile.write( info + "\n" + header )
	outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = -1
	binSize = -1
	outID = ''
	startInd = 0
	percentile = PERCENTILE
	chrmList = None
	
	for i in range( min(6, len(argv)-2) ):
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
		elif argv[i].startswith( '-c=' ):
			chrmList = argv[i][3:].split( ',' )
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
	# end
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	
	fastaIndexStr = argv[startInd]
	
	bedFileStrAr = []
	for j in range( startInd + 1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( fastaIndexStr, bedFileStrAr, numBins, binSize, chrmList, numProc, outID, percentile )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: chrom_plot_bed_combine_pe.py [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-t=percentile] <fasta_index> <bed_file>")
	else:
		parseInputs( sys.argv[1:] )
