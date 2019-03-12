import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *
import bth_util

# Usage: chrom_plot_bed_mid_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] [-t=percentile] <fasta_index> <bed_file> [bed_file]*

NUMPROC=1
NUMBINS=100
PERCENTILE=0.95

def processInputs( bedFileStrAr, fastaIndexStr, labels, outID, chrmList, numProc, numBins, binSize, percentile ):

	if labels == None:
		labels = getSampleNames( bedFileStrAr )
	print( 'Reading FASTA index' )
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	
	#print( chrmDict )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	
	outFileStr = 'chrm_bed'
	if outID != '':
		outFileStr += '_' + outID
	if numBins != -1:
		outFileStr += '_n{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:s}'.format( bth_util.binSizeToStr( binSize ) )
	outFileStr += '.tsv'
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processBEDFile, args=(f, chrmDict, binSize, aNumBins, percentile ) ) for f in bedFileStrAr ]
	bedDictAr = [ p.get() for p in results ]
	
	info = "#from_script:chrom_plot_bed_mid_pe.py; "
	if binSize == -1:
		info += "num_bins:{:d}; ".format( aNumBins )
	else:
		info += "bin_size:{:s}; ".format( bth_util.binSizeToStr( binSize ) )
	info += "num_chrms:{:d}; ".format( len( chrmDict.keys() ) )
	info += "percentile:{:.1f}; ".format( percentile*100 )
	info += "unit:million reads per bin normalized by library size".format( )
	
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, bedDictAr, labels, info )
	
def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
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

def processBEDFile( bedFileStr, chrmDict, binSize, numBins, percentile ):
	print( 'Reading {:s}'.format( os.path.basename(bedFileStr) ) )
	#gffFile = open( gffFileStr, 'r' )
	bedFile = FileBED( bedFileStr )
	bedDict, readCounts = bedFile.getBedDict( middle=True )
	print( 'Processing {:s}'.format( os.path.basename(bedFileStr) ) )
	outDict = processBEDFileHelper( bedDict, readCounts, chrmDict, binSize, numBins )
	# percentile
	thresh = determinePercentile( outDict, percentile )
	outDict['##threshold##'] = thresh
	return outDict

def processBEDFileHelper( bedDict, readCounts, chrmDict, binSize, numBins ):
	outDict = {}
	
	for chrm in chrmDict.keys():
		if binSize == -1:
			outDict[chrm] = [0] * numBins
		else:
			n = math.ceil( chrmDict[chrm] / binSize )
			outDict[chrm] = [0]*n + [-1]*(numBins - n)
	
	# loop through chromosomes
	for chrm in chrmDict.keys():
		# get bed read information
		chrmBedDict = bedDict.get( chrm )
		if chrmBedDict == None:
			print( 'WARNING: {:s} not found in BED file'.format( chrm ) )
			outDict[chrm] = -1
			continue
		if binSize == -1:
			tmpAr = [0] * numBins
			binSize = chrmDict[chrm]
		else:	# have binSize
			n = math.ceil( chrmDict[chrm] / binSize )
			tmpAr = [0]*n + [-1]*(numBins - n)
		for pos in sorted(chrmBedDict.keys()):
			bin = int( pos // binSize )
			tmpAr[bin] += chrmBedDict[pos]
		# normalize by library size
		tmpAr = [ (-1 if x == -1 else x * 1000000 / readCounts) for x in tmpAr ]
		outDict[chrm] = tmpAr	
	return outDict
	
def adjustCounts( gffDict, binSize, chrmDict ):
	for chrm in gffDict.keys():
		getAr = gffDict[ chrm ]
		if binSize == -1:
			divSize = chrmDict[ chrm ]
			newAr = [ ( -1 if x == -1 else float(x) / divSize * 10000 ) for x in getAr ]
		else:
			newAr = getAr
		gffDict[chrm] = newAr
	return gffDict

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

def writeOutput( outFileStr, bedDictAr, labels, info ):
	outFile = open( outFileStr, 'w' )
	header = "sample\tchrm\tbin\tvalue\n"
	outFile.write( info + "\n" + header)
	
	# loop through samples
	for i in range( len( bedDictAr ) ):
		thresh = bedDictAr[i]['##threshold##']
		del bedDictAr[i]['##threshold##']
		# loop through chroms
		for chrm in sorted(bedDictAr[i].keys()):
			chrmAr = bedDictAr[i][chrm]
			if chrmAr == -1:
				continue
			# loop through bin
			for j in range(len( chrmAr )):
				if chrmAr[j] == -1:
					outFile.write( "{:s}\t{:s}\t{:d}\tNA\n".format( labels[i], chrm, j ) )
				elif chrmAr[j] > thresh:
					outFile.write( "{:s}\t{:s}\t{:d}\t{:.3f}\n".format( labels[i], chrm, j, thresh ) )
				else:
					outFile.write( "{:s}\t{:s}\t{:d}\t{:.3f}\n".format( labels[i], chrm, j, chrmAr[j] ) )
	
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = -1
	binSize = -1
	labels = None
	percentile = PERCENTILE
	outID = ''
	startInd = 0
	chrmList = None
	
	for i in range( min(7, len(argv)) ):
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
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			eixt()
	# end
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	fastaIndexStr = argv[startInd]
	
	bedFileStrAr = []
	for j in range( startInd+1, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	if labels != None and len(labels) != len(bedFileStrAr):
		print( "ERROR: number of labels doesn't match number of input files" )
		exit()
	
	
	processInputs( bedFileStrAr, fastaIndexStr, labels, outID, chrmList, numProc, numBins, binSize, percentile )

def printHelp():
	outStr = "Usage: python3 chrom_plot_bed_mid_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] [-t=percentile] <fasta_index> <bed_file> [bed_file]*\n\n"
	outStr += 'Required:\n'
	outStr += 'fasta_index\tfasta index file; ends in .fa.fai\n'
	outStr += 'bed_file\tbed file for ChIP\n\n'
	outStr += 'Optional:\n'
	outStr += '-o=out_id\tidentifier for output [default: num bins or bin size]\n'
	outStr += '-p=num_proc\tnumber of processors to use [default: 1]\n'
	outStr += '-b=num_bins\tbreak chromosomes into this many equally sized bins [default: 100]\n'
	outStr += '\tor\n'
	outStr += '-s=bin_size\tbreak chromosomes into bins of this size [default: NA]\n'
	outStr += '-c=chrm_list\tcomma separated list of chromsomes to include\n\t\t[default: chrms in fasta index]\n'
	outStr += '-l=labels\tcomma-separated list of labels for GFF files [default: gff file name]\n'
	outStr += '-t=percentile\tpercentile to correct against\n'
	print( outStr )
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
