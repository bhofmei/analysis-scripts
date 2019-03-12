import sys, math, glob, multiprocessing, subprocess, os

# Usage: chrom_plot_allc_pe.py [-r] [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_identifier] [-c=chrm_list] [-m=meth_type] <fasta_index> <allc_path> <sample_name> [sample_name]*

NUMPROC=2
NUMBINS=100

def processInputs( allcPath, sampleNamesAr, fastaIndexStr, chrmList, numBins, binSize, isReads, outID, methType, numProc ):
	
	chrmDict = readFastaIndex( fastaIndexStr, chrmList )
	print( 'Read FASTA index.' )
	#print( chrmDict )
	aNumBins, chrmDict = determineNumBins( chrmDict, numBins, binSize )
	
	outFileStr = 'chrm_allc'
	if outID != '':
		outFileStr += '_' + outID
	if numBins != -1:
		outFileStr += '_{:d}'.format( numBins )
	elif binSize != -1:
		outFileStr += '_{:d}'.format( binSize )
	outFileStr += '_{:s}.tsv'.format( methType )
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processSample, args=(allcPath, sampleName, chrmDict, numBins, binSize, methType, isReads) ) for sampleName in sampleNamesAr ]
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	outStrs = [ p.get() for p in results ]
	writeOutput( outFileStr, outStrs )
	print( 'Done' )
	
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

def processSample( allcPath, sampleName, chrmDict, numBins, binSize, methType, isReads ):

	outStr = ''
	
	# loop through chromosomes
	for chrm in sorted(chrmDict.keys()):
		allcFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sampleName, chrm ) )
		print( 'Reading {:s}...'.format( allcFileStr ) )
		chrmAr = readAllcFile( allcFileStr, binSize, chrmDict[chrm], numBins, methType )
		chrmStr = allcArToStr( chrmAr, sampleName, chrm, isReads )
		outStr += chrmStr
	return outStr
	

def readAllcFile( allcFileStr, binSize, chrmDictEntry, numBins, methType ):
	# only for one chromosome -> return array
	
	if binSize == -1:
		outAr = [ [0,0] for x in range(numBins) ]
	else:
		n = math.ceil( chrmDictEntry / binSize )
		outAr = [ [0,0] for x in range(n) ] + [[-1,1] for x in range(numBins - n) ]
	
	allcFile = open( allcFileStr, 'r' )
	
	for line in allcFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		chrm = lineAr[0]
		pos = int( lineAr[1] )
		# only wants C's of correct type
		if methType == 'C' or methType == decodeMethType( lineAr[3] ):
			if binSize == -1:
				bin = int( pos // chrmDictEntry )
			else:
				bin = int( pos // binSize )
			outAr[bin][0] += int( lineAr[4] )
			outAr[bin][1] += int( lineAr[5] )
	allcFile.close()
	return outAr

def decodeMethType( mStr ):
	
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return False
	else:
		return 'CHH'

def allcArToStr( allcAr, sampleName, chrm, isReads ):
	
	outStr = ''
	
	if isReads:
		for i in range(len(allcAr)):
			if allcAr[i][0] == '-1':
				outStr += '{:s}\t{:d}\t{:s}\tNA\tNA\n'.format( chrm, i, sampleName )
			else:
				outStr += '{:s}\t{:d}\t{:s}\t{:d}\t{:d}\n'.format( chrm, i, sampleName, allcAr[i][0], allcAr[i][1] )
	else:
		allcValAr = [ (0 if x[1] == 0 else float(x[0])/float(x[1])) for x in allcAr ]
		for i in range(len(allcAr)):
			if allcValAr[i] == -1:
				outStr += '{:s}\t{:d}\t{:s}\tNA\n'.format( chrm, i, sampleName )
			else:
				outStr += '{:s}\t{:d}\t{:s}\t{:.4f}\n'.format( chrm, i, sampleName, allcValAr[i] )
	return outStr

def writeOutput( outFileStr, outStrAr ):
	
	# header would be: chrm, bin, type, (value)
	
	outFile = open( outFileStr, 'w' )
	for tStr in outStrAr:
		outFile.write( tStr )
	outFile.close()	

def parseInputs( argv ):
	outID = ''
	chrmList = None
	numProc = NUMPROC
	numBins = -1
	binSize = -1
	isReads = False
	methType = 'C'
	startInd = 0
	
	for i in range( min(6, len(argv)-3) ):
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
		elif argv[i] == '-r':
			isReads = True
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			methType = argv[i][3:]
			if methType not in [ 'C', 'CG', 'CHG', 'CHH' ]:
				print( 'ERROR: methylation type must be one of: C, CG, CHG, CHH' )
				exit()
			startInd += 1
	
	if numBins == -1 and binSize == -1:
		numBins = NUMBINS
	
	fastaIndexStr = argv[startInd]
	allcPath = argv[startInd+1]
	if os.path.isdir( allcPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory'.format(allcPath) )
			exit()
	
	sampleNamesAr = []
	for j in range( startInd+2, len(argv) ):
		sampleNamesAr += [ argv[j] ]
	
	processInputs( allcPath, sampleNamesAr, fastaIndexStr, chrmList, numBins, binSize, isReads, outID, methType, numProc )
	

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: chrom_plot_allc_pe.py [-r] [-p=num_proc] [-b=num_bins | -s=bin_size] [-o=output_id] [-c=chrm_list] <fasta_index> <allc_path> <sample_name> [sample_name]*\n-r\t\tprint reads methylated/total rather than weighted methylation\n-p=num_proc\tnumber of processors to use\n-b=num_bins\tnumber of bins to break chromosome into(bin size will differ)\n\t\t[default 100]\nor\n-s=bin_size\tsize (bp) to make each bin (number of bins will differ)\n-o=output_id\tstring to use for output file\n-c=chrm_list\tcomma-separated list of chromsomes to use [default to all listed\n\t\tin fasta index]\nfasta_index\tfasta index file or tab-separated file with chromosomes and\n\t\tlengths\nallc_path\tpath to the allc files\nsample_name\tsample name used in the allc file")
	else:
		parseInputs( sys.argv[1:] )
