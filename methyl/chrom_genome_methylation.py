import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 chrom_genome_methylation.py [-f] [-b] [-p=num_proc] [-o=out_prefix] [-c=chromosomes] [-m=meth_types] <allC_path> <sample_name> [sample_name]*
# reports the weighted methylation (in all contexts) of samples for each chromosome and for the entire genome

NUMPROC=1
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def processInputs( allCPath, sampleNamesAr, numProc, chrmList, fastaIndex, methTypes, outPre, sampleFile, binTest ):

	methTypes.sort()
	
	if binTest == False:
		outFileStr = outPre + '.tsv'
	else:
		outFileStr = outPre + '_binomial.tsv'
	if sampleFile:
		print( 'Reading file of samples names...' )
		sampleNamesAr = readSampleFile( sampleNamesAr[0] )
	if fastaIndex != None:
		chrmList = readFastaIndex( fastaIndex )
	print( '\nMethylation types: {:s}\nChromosomes: {:s}\nUsing binomial test: {:s}\nSamples included: {:s}\nOutput written to: {:s}\n'.format( ' '.join(methTypes), ' '.join(chrmList), str(binTest), ' '.join(sampleNamesAr), outFileStr ) )
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processSample, args=(sample, allCPath, chrmList, methTypes, binTest) ) for sample in sampleNamesAr ]
	sampleDicts = [ p.get() for p in results ]
	
	totalDict = {}
	totalChrmList = chrmList + ['all']
	# set up totalDict
	for chrm in totalChrmList:
		totalDict[chrm] = [None] * len(sampleNamesAr)
	print( "Combining sample's dictionaries..." )
	# add samples
	for i in range(len(sampleNamesAr)):
		# loop through chromosomes
		for chrm in sampleDicts[i].keys():
			totalDict[chrm][i] = sampleDicts[i][chrm]
	
	# write output
	info = '#from_script: chrom_genome_methylation.py; meth_types: {:s}; samples: {:s}; binomial: {:s}'.format( ','.join(methTypes), ','.join(sampleNamesAr), str(binTest) )
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, totalDict, sampleNamesAr, methTypes, info )
		

def readSampleFile( fileStr ):
	
	sampleAr = []
	inFile = open( fileStr, 'r' )
	for line in inFile:
		name = line.rstrip()
		sampleAr += [ name ]
	inFile.close()
	return sampleAr

def readFastaIndex( fastaIndexStr ):
	chrmList = []
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		chrmList += [ lineAr[0] ]
	fastaIndex.close()
	return chrmList

def processSample( sampleName, allCPath, chrmList, methTypes, binTest ):
	
	sampleDict = {}
	genomeMeth = [0] * len( methTypes )
	genomeTotal = [0] * len( methTypes )
	
	# loop through chrms
	for chrm in chrmList:
		allCFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, sampleName, chrm ) )
		print( 'Reading {:s}...'.format( allCFileStr ) )
		chrmMeth, chrmTotal = readAllC( allCFileStr, methTypes, binTest )
		# update genome counts
		genomeMeth = addToArray( genomeMeth, chrmMeth )
		genomeTotal = addToArray( genomeTotal, chrmTotal )
		# compute methylation
		sampleDict[ chrm ] = computeMethylation( chrmMeth, chrmTotal )
	# add genome
	sampleDict[ 'all' ] = computeMethylation( genomeMeth, genomeTotal )
	return sampleDict

def readAllC( allCFileStr, methTypes, binTest ):
	
	methC = [0] * len( methTypes )
	totalC = [0] * len( methTypes )
	
	allCFile = open( allCFileStr, 'r' )
	for line in allCFile:
		# header
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		
		# skip unmethylated positions when binTest is true
		if binTest and lineAr[6] == '0':
			continue
		mType = decodeMethType( lineAr[3] )
		if mType in methTypes:
			ind = methTypes.index( mType )
			methC[ind] += int( lineAr[4] )
			totalC[ind] += int( lineAr[5] )
		if 'C' in methTypes:
			ind = methTypes.index( 'C' )
			methC[ind] += int( lineAr[4] )
			totalC[ind] += int( lineAr[5] )
	# end for
	allCFile.close()
	return methC, totalC

def decodeMethType( mStr ):
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'

def addToArray( oldArray, tmpArray ):
	if len(oldArray) != len(tmpArray):
		print( 'ERROR: array lengths differ' )
		return
	for j in range( len(oldArray) ):
		oldArray[j] += tmpArray[j]
	return oldArray

def computeMethylation( methAr, totalAr ):
	
	outAr = []
	for i in range(len(methAr)):
		weighted = float( methAr[i] ) / float( totalAr[i] )
		outAr += [ '{:.4f}'.format( weighted ) ]
	return outAr

def writeOutput( outFileStr, totalDict, sampleNamesAr, methTypes, info ):
	
	outFile = open( outFileStr, 'w' )
	# header
	outStr = info + '\nsample\tchrm\t' + '\t'.join( methTypes ) + '\n'
	outFile.write( outStr )
	
	# loop through chromosomes
	for chrm in sorted(totalDict.keys()):
		tmpMatrix = totalDict[chrm]
		# loop through samples
		for i in range(len(sampleNamesAr)):
			outStr = '{:s}\t{:s}\t{:s}\n'.format( sampleNamesAr[i], chrm, '\t'.join( tmpMatrix[i] ) )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	
	outPre = 'out'
	numProc = NUMPROC
	chrmList = None
	fastaIndex = None
	sampleFile = False
	binTest = False
	methTypes = ['CG','CHG','CHH','C']
	startInd = 0
	
	for i in range( min(6, len(argv)-2) ):
		
		if argv[i] == '-f':
			sampleFile = True
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i].startswith( '-o='):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c='):
			if fastaIndex != None:
				print( 'ERROR: cannot specify chromosome list and fasta index file' )
				exit()
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-cf=' ):
			if chrmList != None:
				print( 'ERROR: cannot specify chromosome list and fasta index file' )
				exit()
			fastaIndex = argv[i][4:]
			startInd += 1
		elif argv[i].startswith( '-m='):
			methTypes = argv[i][3:].split(',')
			startInd += 1
		elif argv[i] == '-b':
			binTest = True
			startInd += 1
	# end for
	if chrmList == None and fastaIndex == None:
		chrmList = CHRMLIST
	
	allCPath = argv[startInd]
	if os.path.isdir( allCPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory'.format(allCPath) )
			exit()
			
	sampleNamesAr = []
	for j in range( startInd+1, len(argv) ):
		sampleNamesAr += [ argv[j] ]
		if j == 1 and sampleFile:
			print( 'ERROR: only specify one file with sample names' )
			exit()
	processInputs( allCPath, sampleNamesAr, numProc, chrmList, fastaIndex, methTypes, outPre, sampleFile, binTest )
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 chrom_genome_methylation.py [-f] [-b] [-p=num_proc] [-o=out_prefix] [-c=chromosomes | -cf=fasta_index] [-m=meth_types] <allC_path> <sample_name> [sample_name]")
		print("-f\tsamples are listed in file (1 line per sample)\n-c\tcomma-separated list of chromosomes to analyze [default Chr1-5]\nor\n-cf=fasta_index\tuse chromosomes in fasta index file\n-m\tcomma-separated list methylation contexts to analyze [default CG,CHG,CHH,C]\n-p\tnumber of processors to use [default 2]\n-o\tprefix to use for output file[default 'out']\n-b\tcompute weighted methylation for positions passing binomial test only")
	else:
		parseInputs( sys.argv[1:] )
