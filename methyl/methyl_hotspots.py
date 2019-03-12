import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: methyl_hotspots.py [-h] [-o=out_prefix] [-b=bin_size] [-ex=exclude_names] [-in=include_names] <in_file>

BINSIZE=1000000
OUTPRE='out'
MINSIZE=1000000

def processInputs(inFileStr, binSize, outPre, exList, inList, minSize ):
	
	# read file
	coordDict, maxCoord = readFile( inFileStr, minSize )
	# clean dictionary and get sample names
	if exList != None or inList != None:
		coordDict = cleanDict( coordDict, exList, inList )
	sampleAr = sorted(list( coordDict.keys() ))
	# count
	countAr = countRecomb( coordDict, binSize, maxCoord )
	#print(len(countAr))
	# header
	info = "binsize:{:d}, samples:{:s}".format( binSize, ','.join(sampleAr) )
	
	# write output
	outFileStr = outPre + "_hotspots.tsv"
	writeOutput( outFileStr, countAr, binSize, info )

def readFile( inFileStr, minSize ):
	# getting coordinate values for samples
	maxCoord = 0
	coordDict = {}
	inFile = open( inFileStr, 'r' )
	for line in inFile:
		if line.startswith( "[" ):
			continue
		lineAr = line.rstrip().split( ' ' )
		samp = lineAr[0]
		coordDict[samp] = []
		prevVal = None
		for j in range( 1, len(lineAr), 3 ):
			val = int(lineAr[j])
			pos1 = int(lineAr[j+1])
			pos2 = int(lineAr[j+2])
			if pos2>maxCoord:
				maxCoord = pos2
			if pos2 - pos1 + 1 >= minSize:
				if val == prevVal:
					coordDict[samp][-1] = pos2
				else:
					coordDict[samp] += [ pos2 ]
				prevVal = val
	inFile.close()
	return coordDict, maxCoord

def cleanDict( coordDict, exList, inList ):
	sampleNames = list( coordDict.keys() )
	newDict = {}
	# include list overrides exclude list
	if inList != None:
		for sample in inList:
			if sample in sampleNames:
				newDict[sample] = coordDict[sample]
	else:
		for sample in sampleNames:
			if sample not in exList:
				newDict[sample] = coordDict[sample]
	return newDict

def countRecomb( coordDict, binSize, maxCoord ):
	
	outAr = [0] * int(math.ceil(maxCoord/float(binSize)))
	# loop through samples
	for key in coordDict.keys():
		# loop though all but last coordinate
		keyAr = coordDict[key]
		for i in range(len(keyAr)-1):
			bin = ( keyAr[i]-1 ) // binSize
			outAr[bin] += 1
	return outAr

def writeOutput( outFileStr, countAr, binSize, info ):
	header = "#bin\tcount\n#" + info + "\n"
	outFile = open( outFileStr, 'w' )
	outFile.write( header )
	
	for i in range(len(countAr)):
		outFile.write( '{:d}\t{:d}\n'.format( i * binSize + 1, countAr[i] ) )
	outFile.close()

def printHelp():
	out = "Usage: methyl_hotspots.py [-h] [-o=out_prefix] [-b=bin_size] [-ex=exclude_names] [-in=include_names] <in_file>"
	out += "Required:\nin_file\tinput file unformatted\n"
	out += "Optional:\n"
	out += "-h\tprint this help screen and exit\n"
	out += "-o=out_prefix\tprefix to use for output file\n\t[default:'out']\n"
	out += "-b=bin_size\tbin size to use for counting\n\t[default:100KB]\n"
	out += "-ex=exclude_names\tcomma-separated list of sample names to exclude\n"
	out += "-in=include_names\tcomma-separated list of the sample names to only include\n"
	out += "if -in and -ex both set, uses include list\n"
	return out
	
def parseInputs( argv ):
	binSize = BINSIZE
	outPre = OUTPRE
	minSize = MINSIZE
	exList = None
	inList = None
	startInd = 0
	
	for i in range(min(6,len(argv)-1)):
		if argv[i] == "-h":
			print(printHelp())
			exit()
		elif argv[i].startswith('-b='):
			try:
				binStr = argv[i][3:].lower()
				if binStr.endswith( 'k' ):
					binSize = int( binStr[:-1] ) * 1000
				elif  binStr.endswith( 'kb' ):
					binSize = int( binStr[:-2] ) * 1000
				elif binStr.endswith( 'm' ):
					binSize = int( binStr[:-1] ) * 1000000
				elif binStr.endswith( 'mb' ):
					binSize = int( binStr[:-2] ) * 1000000
				else:
					binSize = int( binStr )
				startInd += 1
			except ValueError:
				print( 'ERROR: bin size must be an integer' )
				exit()
		elif argv[i].startswith('-o='):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-ex='):
			exList = argv[i][4:].split(',')
			startInd += 1
		elif argv[i].startswith('-in='):
			inList = argv[i][4:].split(',')
			startInd += 1
		elif argv[i].startswith('-m='):
			try:
				minStr = argv[i][3:].lower()
				if minStr.endswith( 'k' ):
					minSize = int( minStr[:-1] ) * 1000
				elif minStr.endswith( 'kb' ):
					minSize = int( minStr[:-2] ) * 1000
				elif minStr.endswith( 'm' ):
					minSize = int( minStr[:-1] ) * 1000000
				elif minStr.endswith( 'mb' ):
					minSize = int( minStr[:-2] ) * 1000000
				else:
					minSize = int( minStr )
				startInd += 1
			except ValueError:
				print( 'ERROR: min size must be an integer' )
				exit()
			
	inFileStr = argv[startInd]
	processInputs( inFileStr, binSize, outPre, exList, inList, minSize )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print (printHelp())
	else:
		parseInputs( sys.argv[1:] )
