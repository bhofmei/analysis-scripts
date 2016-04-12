import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: chrom_plot_gff_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] <fasta_index> <gff_file> [gff_file]*

NUMPROC=2
NUMBINS=100

def processInputs(gffFileStrAr, fastaIndexStr, labels, outID, chrmList, numProc, numBins, binSize):

	if labels == None:
		labels = getSampleNames( gffFileStrAr )
		
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	print( 'Read FASTA index.' )
	#print( chrmDict )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	
	outFileStr = 'chrm_gff'
	if outID != '':
		outFileStr += '_' + outID
	if numBins != -1:
		outFileStr += '_n{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:d}'.format( binSize )
	outFileStr += '.tsv'
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processGFFFile, args=(f, chrmDict, binSize, aNumBins ) ) for f in gffFileStrAr ]
	gffDictAr = [ p.get() for p in results ]
	
	info = "#from_script:chrom_plot_gff_pe.py;"
	if binSize == -1:
		info += "num_bins:{:d}".format( aNumBins )
	else:
		info += "bin_size:{:d}".format( binSize )
	info += "; num_chrms:{:d}".format( len( chrmDict.keys() ) )
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, gffDictAr, labels, info )
	
def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	leftIndex = fileStr.rfind('/')
	rightIndex = fileStr.rfind('.')
	sampleName = fileStr[leftIndex+1:rightIndex]
	sampleNamesAr += [ sampleName ]
	return sampleNamesAr

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


def processGFFFile( gffFileStr, chrmDict, binSize, numBins ):
	gffDict = {}
	
	for chrm in chrmDict.keys():
		if binSize == -1:
			gffDict[chrm] = [0] * numBins
		else:
			n = math.ceil( chrmDict[chrm] / binSize )
			gffDict[chrm] = [0]*n + [-1]*(numBins - n)
	print( 'Reading {:s}...'.format( gffFileStr ) )
	gffFile = open( gffFileStr, 'r' )
	
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		chrm = lineAr[0]
		success = chrmDict.get( chrm )
		featType = lineAr[2]
		if success == None or featType not in ["gene", "similarity"]:
			continue
		if binSize == -1:
			bStart = int( start // chrmDict[chrm] )
			bEnd = int( end // chrmDict[ chrm] )
		else:
			bStart = int( start // binSize )
			bEnd = int( end // binSize )
		#print( start, end, bStart, bEnd )
		for bin in range( bStart, bEnd + 1):
			gffDict[chrm][bin] += 1
	
	gffFile.close()
	# adjust by bin size
	#gffDict = adjustCounts( gffDict, binSize, chrmDict )
	return gffDict
	
def adjustCounts( gffDict, binSize, chrmDict ):
	for chrm in gffDict.keys():
		getAr = gffDict[ chrm ]
		if binSize == -1:
			divSize = chrmDict[ chrm ]
		else:
			divSize = binSize
		newAr = [ ( -1 if x == -1 else float(x) / divSize * 100000 ) for x in getAr ]
		gffDict[chrm] = newAr
	return gffDict

def writeOutput( outFileStr, gffDictAr, labels, info ):
	outFile = open( outFileStr, 'w' )
	header = "#sample\tchrm\tbin\tvalue\n"
	outFile.write( info + "\n" + header)
	
	# loop through samples
	for i in range( len( gffDictAr ) ):
		# loop through chroms
		for chrm in sorted(gffDictAr[i].keys()):
			chrmAr = gffDictAr[i][chrm]
			# loop through bin
			for j in range(len( chrmAr )):
				if chrmAr[j] == -1:
					outFile.write( "{:s}\t{:s}\t{:d}\tNA\n".format( labels[i], chrm, j ) )
				else:
					outFile.write( "{:s}\t{:s}\t{:d}\t{:.3f}\n".format( labels[i], chrm, j, chrmAr[j] ) )
	
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = -1
	binSize = -1
	labels = None
	outID = ''
	startInd = 0
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
		elif argv[i].startswith( '-l=' ):
			labels = argv[i][3:].split( ',' )
			startInd += 1
	# end
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	fastaIndexStr = argv[startInd]
	
	gffFileStrAr = []
	for j in range( startInd+1, len(argv) ):
		gffFileStrAr += [ argv[j] ]
	
	if len(labels) != len( gffFileStrAr):
		print( "ERROR: number of labels doesn't match number of input files" )
		exit
	
	processInputs( gffFileStrAr, fastaIndexStr, labels, outID, chrmList, numProc, numBins, binSize )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: chrom_plot_gff_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] <fasta_index> <gff_file> [gff_file]*")
	else:
		parseInputs( sys.argv[1:] )
