import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 adjust_meta_intergenic.py [-e] [-o=out_file] [-b=num_bins[,num_bins_stream]] <-s | -d>  <metaplot_file>

def processInputs( metaFileStr, outFileStr, adMode, numBins, numBinsStream, expFile ):
	
	if metaFileStr.endswith( '.tsv' ):
		isTSV = True
	else:
		isTSV = False
	if outFileStr == None:
		outFileStr = getOutFileName( metaFileStr, adMode, isTSV )
	if numBins == -1:
		numBins = getNumBins( metaFileStr )
	print( 'Reading {:s}...'.format( metaFileStr ) )
	if expFile:
		metaDict, infoLine = readMetaFileExpression( metaFileStr, numBins, numBinsStream, isTSV )
	else:
		metaDict, infoLine = readMetaFile( metaFileStr, numBins, numBinsStream, isTSV )
	print( 'Adjusting values using {:s} mode...'.format( adMode ) )
	adjustDict = adjustMetaDict( metaDict, adMode )
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, adjustDict, infoLine, expFile, isTSV )
	

def getOutFileName( metaFileStr, adMode, isTSV ):
	
	rInd = metaFileStr.rfind( '.' )
	outFileStr = metaFileStr[:rInd] + '_adjust'
	outFileStr += ('-s' if adMode == 'subtract' else '-d' )
	outFileStr += ( '.tsv' if isTSV else '.csv' )
	return outFileStr

def getNumBins( metaFileStr ):
	'''
		to be used when meta file is output of metaplot_from_bed_allc_pe.py
	'''
	
	fileName = os.path.basename( metaFileStr )
	# remove the leading 'meta_'
	fileName = fileName[5:]
	lInd = fileName.find( '_' )
	try:
		numBins = int( fileName[:lInd] )
	except ValueError:
		print( 'ERROR: specify number of bins and optionally number of bins for stream' )
		exit()
	return numBins


def readMetaFile( metaFileStr, numBins, numBinsStream, isTSV ):
	
	# set up a dictionary
	# {sample: [ [intergenic],[genic] ], sample: [ [intergenic], [genic] ] }
	# where intergenic and genic are arrays of meta values
	
	metaFile = open( metaFileStr, 'r' )
	
	storeArray = []
	curSample = None
	metaDict = {}
	infoLine = ''
	delim = ( '\t' if isTSV else ',' )
	
	for line in metaFile:
		if line.startswith( 'bin' ):
			continue
		elif line.startswith( '#' ):
			infoLine += line
			continue
		lineAr = line.rstrip().split( delim )
		# (0) bin number (1) value (2) sample
		# first sample
		if curSample == None:
			curSample = lineAr[2]
		# change sample -> save to dictionary
		elif curSample != lineAr[2]:
			# see if need to get num bins stream
			if numBinsStream == -1:
				numBinsStream = getNumBinsStream( storeArray, numBins )
			# save to dictionary
			metaDict[ curSample ] = splitArray( storeArray, numBins, numBinsStream )
			storeArray = []
			curSample = lineAr[2]
		# save bin
		storeArray += [ float( lineAr[1] ) ]
	
	metaDict[curSample] = splitArray( storeArray, numBins, numBinsStream )
	metaFile.close()
	
	return metaDict

def readMetaFileExpression( metaFileStr, numBins, numBinsStream, isTSV ):
	
	# set up a dictionary
	# {sample: [ [intergenic],[genic] ], sample: [ [intergenic], [genic] ] }
	# where intergenic and genic are arrays of meta values
	
	metaFile = open( metaFileStr, 'r' )
	
	storeArray = []
	curSample = None
	metaDict = {}
	infoLine = ''
	delim = ( '\t' if isTSV else ',' )
	
	for line in metaFile:
		if line.startswith( 'sample' ):
			continue
		elif line.startswith( '#' ):
			infoLine += line
			continue
		lineAr = line.rstrip().split( delim )
		s = lineAr[0]+delim+lineAr[1]
		# (0) sample (1) expression (2) bin (3) value
		# first sample
		if curSample == None:
			curSample = s
		# change sample -> save to dictionary
		elif curSample != s:
			# see if need to get num bins stream
			if numBinsStream == -1:
				numBinsStream = getNumBinsStream( storeArray, numBins )
			# save to dictionary
			metaDict[ curSample ] = splitArray( storeArray, numBins, numBinsStream )
			storeArray = []
			curSample = s
		# save bin
		storeArray += [ float( lineAr[3] ) ]
	
	metaDict[curSample] = splitArray( storeArray, numBins, numBinsStream )
	metaFile.close()
	
	return metaDict, infoLine

def getNumBinsStream( storeArray, numBins ):
	return int ((len( storeArray ) - numBins) / 2)
	
def splitArray( inArray, numBins, numBinsStream ):
	
	interAr = inArray[:numBinsStream] + inArray[(numBinsStream+numBins):]
	geneAr = inArray[numBinsStream:(numBinsStream + numBins)]
	return [ interAr, geneAr ]

def adjustMetaDict( metaDict, adMode ):
	
	outDict = {}
	for key in metaDict.keys():
		#print (key)
		interAr = metaDict[key][0]
		geneAr = metaDict[key][1]
		if adMode == 'subtract':
			adAr = adjustSubtract( interAr, geneAr )
		elif adMode == 'divide':
			adAr = adjustDivide( interAr, geneAr )
		outDict[ key ] = adAr
	return outDict

def adjustSubtract( interAr, geneAr ):
	
	interAve = math.fsum( interAr ) / len( interAr )
	print( interAve )
	
	interArAd = [ x - interAve for x in interAr ]
	geneArAd = [ x - interAve for x in geneAr ]
	
	interInd = int( len(interArAd) / 2 )
	
	outAr = interArAd[:interInd] + geneArAd + interArAd[interInd:]
	return outAr

def adjustDivide( interAr, geneAr ):
	
	interAve = math.fsum( interAr ) / len( interAr )
	print( interAve )
	
	interArAd = [ x / interAve for x in interAr ]
	geneArAd = [ x / interAve for x in geneAr ]
	
	interInd = int( len(interArAd) / 2 )
	
	outAr = interArAd[:interInd] + geneArAd + interArAd[interInd:]
	return outAr

def writeOutput( outFileStr, adjustDict, infoLine, expFile, isTSV ):
	
	outFile = open( outFileStr, 'w' )
	# header
	delim = ( '\t' if isTSV else ',' )
	if expFile:
		outFile.write( 'sample{:s}expression{:s}bin{:s}value\n'.format( delim, delim, delim ))
	else:
		outFile.write( 'bin{:s}count{:s}sample\n'.format( delim, delim ) )
	if infoLine != '':
		outFile.write( infoLine )
	
	for key in sorted(adjustDict.keys()):
		adAr = adjustDict[key]
		for i in range(len(adAr)):
			if expFile:
				outFile.write( '{:s}{:s}{:d}{:s}{:f}\n'.format( key, delim, i+1, delim, adAr[i] ) )
			else:
				outFile.write( '{:d}{:s}{:f}{:s}{:s}\n'.format( i+1, delim, adAr[i], delim, key ) )
	
	outFile.close()

def parseInputs( argv ):
	# Usage: python3.4 adjust_meta_intergenic.py [-e] [-o=out_file] [-b=num_bins[,num_bins_stream]] <-s | -d>  <metaplot_file>
	
	outFileStr = None
	adMode = None
	numBins = -1
	numBinsStream = -1
	expFile = False
	startInd = 0
	
	for i in range( min(4,len(argv)) ):
		if argv[i].startswith( '-o=' ):
			outFileStr = argv[i][3:]
		elif argv[i] == '-e':
			expFile = True
			startInd += 1
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
		# handle adjustment modes
		elif argv[i] == '-s' and adMode == None:
			adMode = 'subtract'
		elif argv[i] == '-d' and adMode == None:
			adMode = 'divide'
		elif argv[i] == '-s' or argv[i] == '-d':
			print( 'ERROR: specify only one adjustment mode type' )
			exit()
	metaFileStr = argv[-1]
	
	processInputs( metaFileStr, outFileStr, adMode, numBins, numBinsStream, expFile )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 adjust_meta_intergenic.py [-e] [-o=out_file] [-b=num_bins[,num_bins_stream]] <-s | -d>  <metaplot_file>")
	else:
		parseInputs( sys.argv[1:] )
