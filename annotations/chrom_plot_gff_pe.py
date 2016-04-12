import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import bth_util

# Usage: chrom_plot_gff_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] [-t=calc_type] <fasta_index> <gff_file> [gff_file]*

NUMPROC=1
NUMBINS=100

def processInputs(gffFileStrAr, fastaIndexStr, labels, calcType, outID, chrmList, numProc, numBins, binSize):

	if labels == None:
		labels = getSampleNames( gffFileStrAr )
		
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	print( 'Read FASTA index.' )
	#print( chrmDict )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	
	outFileStr = 'chrm_gff'
	if outID != '':
		outFileStr += '_' + outID
	outFileStr += '_' + ( 'length' if calcType == 'l' else 'number' )
	if numBins != -1:
		outFileStr += '_n{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:s}'.format( bth_util.binSizeToStr( binSize ) )
	outFileStr += '.tsv'
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processGFFFile, args=(f, chrmDict, binSize, aNumBins, calcType ) ) for f in gffFileStrAr ]
	gffDictAr = [ p.get() for p in results ]
	
	info = "#from_script:chrom_plot_gff_pe.py; "
	if binSize == -1:
		info += "num_bins:{:d}".format( aNumBins )
	else:
		info += "bin_size:{:s}".format( bth_util.binSizeToStr( binSize ) )
	info += "; num_chrms:{:d};".format( len( chrmDict.keys() ) )
	info += " unit:{:s} per ".format( 'number' if calcType == 'n' else 'total bp' )
	if binSize == -1:
		info += '10kb'
	elif (binSize // 1000000) > 0:
		info += str( binSize / 1000000 ) + 'mbp'
	elif (binSize // 1000) > 0:
		info += str( binSize / 1000 ) + 'kbp'
	else:
		info += str( binSize ) + 'bp'
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, gffDictAr, labels, info )
	
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

def processGFFFile( gffFileStr, chrmDict, binSize, numBins, calcType ):
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
		chrm = lineAr[0]
		start = int(lineAr[3])
		end = int(lineAr[4])
		success = chrmDict.get( chrm )
		featType = lineAr[2]
		if success == None or featType not in ["gene", "similarity"]:
			continue
		if binSize == -1:
			bStart = int( start // chrmDict[chrm] )
			bEnd = int( end // chrmDict[chrm] )
		else:
			bStart = int( start // binSize )
			bEnd = int( end // binSize )
		#print( start, end, bStart, bEnd )
		if calcType == 'n':
			for bin in range( bStart, bEnd + 1):
				gffDict[chrm][bin] += 1
		else:
			# one bin
			if bStart == bEnd:
				#print( '1.', chrm, start, end, 'a)', bStart, end - start + 1 )
				gffDict[chrm][bStart] += end - start + 1
			# two bins
			elif bStart + 1 == bEnd:
				#print( '2.',chrm, start, end, 'a)',bStart, ( bEnd * binSize ) - start + 1, 'b)', bEnd, end - ( bEnd * binSize ) + 1 )
				gffDict[chrm][bStart] += ( bEnd * binSize ) - start + 1
				gffDict[chrm][bEnd] += end - ( bEnd * binSize ) + 1
			# 3+ bins
			else:
				#print( '3.',chrm, start, end, 'a)',bStart ( (bStart+1) * binSize ) - start + 1, 'b)', bStart+1,'-',bEnd-1,binSize, 'c)',bEnd, end - ( bEnd * binSize ) + 1 )
				gffDict[chrm][bStart] += ( (bStart+1) * binSize ) - start + 1
				for bin in range( bStart+1, bEnd ):
					gffDict[chrm][bin] += binSize
				gffDict[chrm][bEnd] += end - ( bEnd * binSize ) + 1
	# end for line
	
	gffFile.close()
	# adjust by bin size
	#if featType != gene
	gffDict = adjustCounts( gffDict, binSize, chrmDict )
	return gffDict
	
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

def writeOutput( outFileStr, gffDictAr, labels, info ):
	outFile = open( outFileStr, 'w' )
	header = "sample\tchrm\tbin\tvalue\n"
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
	calcType = 'n'
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
			x = argv[i][3:]
			if x.lower() not in ['l','n']:
				print( 'ERROR: {:s} is not a valid calculation type'.format(x) )
				exit()
			calcType = x.lower()
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			eixt()
	# end
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	fastaIndexStr = argv[startInd]
	
	gffFileStrAr = []
	for j in range( startInd+1, len(argv) ):
		gffFileStrAr += [ argv[j] ]
	
	if labels != None and (len(labels) != len(gffFileStrAr)):
		print( "ERROR: number of labels doesn't match number of input files" )
		exit()
	
	
	processInputs( gffFileStrAr, fastaIndexStr, labels, calcType, outID, chrmList, numProc, numBins, binSize )

def printHelp():
	outStr = "Usage: python3 chrom_plot_gff_pe.py [-o=out_id] [-p=num_proc] [-b=num_bins | -s=bin_size] [-c=chrm_list] [-l=labels] [-t=calc_type] <fasta_index> <gff_file> [gff_file]*\n\n"
	outStr += 'Required:\n'
	outStr += 'fasta_index\tfasta index file; ends in .fa.fai\n'
	outStr += 'gff_file\tgff file of features\n\n'
	outStr += 'Optional:\n'
	outStr += '-o=out_id\tidentifier for output [default: num bins or bin size]\n'
	outStr += '-p=num_proc\tnumber of processors to use [default: 1]\n'
	outStr += '-b=num_bins\tbreak chromosomes into this many equally sized bins [default: 100]\n'
	outStr += '\tor\n'
	outStr += '-s=bin_size\tbreak chromosomes into bins of this size [default: NA]\n'
	outStr += '-c=chrm_list\tcomma separated list of chromsomes to include\n\t\t[default: chrms in fasta index]\n'
	outStr += '-l=labels\tcomma-separated list of labels for GFF files [default: gff file name]\n'
	outStr += '-t=calc_type\tsingle letter indicating calculation type;\n\t\tn - total number of features per bin [default]\n\t\tl - total length of features per bin\n'
	print( outStr )
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
