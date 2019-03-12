import sys, math, glob, multiprocessing, subprocess, os


# Usage: python3.4 chrom_plot_bed_pe.py [-peaks] [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-t=percentile] <fasta_index> <bed_file> [bed_file]*
# splits each chromosome into bins and reports the frequency of chip reads within the bins; each chromosome can have same number of bins or bins can be same number of base pairs

NUMPROC=2
NUMBINS=100
PERCENTILE=0.95

def processInputs( fastaIndexStr, bedFileStrAr, numBins, binSize, chrmList, numProc, outID, percentile, callPeaks ):
	
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	print( 'Read FASTA index.' )
	#print( chrmDict )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	#print( aNumBins )
	#print( chrmDict )
	
	outFileStr = 'chrm_bed'
	peakFileStr = 'peaks'
	if outID != '':
		outFileStr += '_' + outID
		peakFileStr += '_' + outID
	if numBins != -1:
		outFileStr += '_{:d}'.format( numBins )
		peakFileStr += '_{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:d}'.format( binSize )
		peakFileStr += '_{:d}'.format( binSize )
	outFileStr += '_{:d}'.format( int( percentile * 100 ) )
	peakFileStr += '_{:d}'.format( int( percentile * 100 ) )
	outFileStr += '.tsv'
	peakFileStr += '.tsv'
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processBedFile, args=(f, chrmDict, aNumBins, binSize, percentile) ) for f in bedFileStrAr ]
	sStrs = [ p.get() for p in results ]
	outStrs = [ x[0] for x in sStrs ]
	peakStrs = [ x[1] for x in sStrs ]
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, outStrs )
	
	if callPeaks:
		print( 'Writing peaks to {:s}...'.format( peakFileStr ) )
		writeOutput( peakFileStr, peakStrs )
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

def processBedFile( bedFileStr, chrmDict, numBins, binSize, percentile ):
	
	sName = os.path.basename( bedFileStr )
	ind = sName.rfind( '.' )
	sampleName = sName[:ind]
	
	print( 'Processing {:s}...'.format( sampleName ) )
	bedDict, countReads = readBedFile( bedFileStr, chrmDict, numBins, binSize )
	thresh = determinePercentile( bedDict, percentile )
	bedStr, peaksStr = bedDictToStr( bedDict, countReads, sampleName, thresh )
	
	return (bedStr, peaksStr)
		
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
	return bedDict, countReads

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

def bedDictToStr( bedDict, countReads, name, thresh ):
	
	outStr = ''
	peaksStr = ''
	# loop through chroms
	for chrm in sorted( bedDict.keys() ):
		tmpAr = bedDict[chrm]
		# loop through bins
		for i in range( len(tmpAr) ):
			if tmpAr[i] == -1:
				outStr += '{:s}\t{:d}\t{:s}\tNA\n'.format( chrm, i, name )
			elif tmpAr[i] > thresh:
				outStr += '{:s}\t{:d}\t{:s}\t{:.4f}\n'.format( chrm, i, name, (float(thresh) / countReads * 1000) )
				peaksStr += '{:s}\t{:d}\t{:s}\t{:.4f}\n'.format( chrm, i, name, (float(tmpAr[i]) / countReads * 1000) )
			else:
				outStr += '{:s}\t{:d}\t{:s}\t{:.4f}\n'.format( chrm, i, name, (float(tmpAr[i]) / countReads * 1000) )
	return outStr, peaksStr

def writeOutput( outFileStr, outStrAr ):
	
	# header would be: chrm, bin, type, (value)
	
	outFile = open( outFileStr, 'w' )
	for tStr in outStrAr:
		outFile.write( tStr )
	outFile.close()
	

def parseInputs( argv ):
	# Usage: python3.4 chrom_plot_bed_pe.py [-peaks] [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-t=percentile] <fasta_index> <bed_file> [bed_file]*
	
	numProc = NUMPROC
	numBins = -1
	binSize = -1
	outID = ''
	callPeaks = False
	startInd = 0
	percentile = PERCENTILE
	chrmList = None
	
	for i in range( min(7, len(argv)-2) ):
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
				if binStr.endswith( 'k' ):
					binSize = int( binStr[:-1] ) * 1000
				elif binStr.endswith( 'm' ):
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
		elif argv[i] == '-peaks':
			callPeaks = True
			startInd += 1
	# end
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	
	fastaIndexStr = argv[startInd]
	
	bedFileStrAr = []
	for j in range( startInd + 1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( fastaIndexStr, bedFileStrAr, numBins, binSize, chrmList, numProc, outID, percentile, callPeaks )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 chrom_plot_bed_pe.py [-peaks] [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-t=percentile] <fasta_index> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
