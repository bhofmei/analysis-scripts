import sys, math, glob, multiprocessing, subprocess, os, gzip

# sage: python3.4 genome_methylation.py [-q] [-f] [-b] [-p=num_proc] [-o=out_prefix] [-i=include_chrms | -x=exclude_chrms] [-m=meth_types] <allC_path> <sample_name> [sample_name]*
# reports the weighted methylation (in all contexts) of samples genome-wide

NUMPROC=1

def processInputs( allCPath, sampleNamesAr, numProc, includeList, excludeList, methTypes, outPre, sampleFile, binTest, isPrint ):

	methTypes.sort()
	
	if binTest == False:
		outFileStr = outPre + '.tsv'
	else:
		outFileStr = outPre + '_binomial.tsv'
	if sampleFile:
		sampleNamesAr = readSampleFile( sampleNamesAr[0] )
	if isPrint:
		print( 'Methylation types: {:s}\nUsing binomial test: {:s}\nIncludes chromsomes: {:s}\nExclude chromsomes: {:s}\nSamples included: {:s}\nOutput written to: {:s}\n'.format( ' '.join(methTypes), str(binTest), ' '.join(includeList),' '.join(excludeList),' '.join(sampleNamesAr), outFileStr ) )
	
	# check for all allC files before doing any work
	sampleFilesAr = []
	for sample in sampleNamesAr:
		checked  = checkFiles( allCPath, sample )
		if checked == False:
			exit()
		else:
			sampleFilesAr += [checked]
	
	if isPrint:
		print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processSample, args=(sampleFile, includeList, excludeList, methTypes, binTest, isPrint) ) for sampleFile in sampleFilesAr ]
	sampleMethylAr = [ p.get() for p in results ]
	
	# write output
	info = '#from_script: genome_methylation.py; meth_types: {:s}; samples: {:s}; binomial: {:s}'.format( ','.join(methTypes), ','.join(sampleNamesAr), str(binTest) )
	if isPrint:
		print( 'Writing output to {:s}'.format( outFileStr ) )
	writeOutput( outFileStr, sampleMethylAr, sampleNamesAr, methTypes, info )
		

def readSampleFile( fileStr ):
	
	sampleAr = []
	inFile = open( fileStr, 'r' )
	for line in inFile:
		name = line.rstrip()
		sampleAr += [ name ]
	inFile.close()
	return sampleAr

def checkFiles( allcPath, sample ):
	filePath = os.path.normpath('{:s}/allc_{:s}.tsv'.format( allcPath, sample ) )
	if os.path.isfile( filePath ):
		return filePath
	elif os.path.isfile( filePath+'.gz' ):
		return filePath + '.gz'
	else:
		print( 'ERROR: allC file for sample {:s} not found, {:s}'.format( sample, filePath ) )
		return False

def processSample( allCFileStr, includeList, excludeList, methTypes, binTest, isPrint ):
	
	methC = [0] * len( methTypes )
	totalC = [0] * len( methTypes )
	if isPrint:
		print('Reading {:s}'.format(os.path.basename(allCFileStr)))
	if allCFileStr.endswith('.gz'):
		allCFile = gzip.open(allCFileStr, 'rt' )
	else:
		allCFile = open( allCFileStr, 'r' )
	
	for line in allCFile:
		# header
		if line.startswith( 'chr' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) < 7:
			#print(lineAr)
			continue
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		chrm = lineAr[0]
		
		# skip unmethylated positions when binTest is true
		if binTest and lineAr[6] == '0':
			continue
		# skip chrms in exclude list
		if chrm in excludeList:
			continue
		# if include specified and not in it, skip
		if len(includeList) > 0 and chrm in includeList:
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
	
	# compute weighted and return
	return computeMethylation(methC, totalC)

def decodeMethType( mStr ):
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'

def computeMethylation( methAr, totalAr ):
	
	outAr = []
	for i in range(len(methAr)):
		weighted = float( methAr[i] ) / float( totalAr[i] )
		outAr += [str(methAr[i])] + [str(totalAr[i])] + [ '{:.4f}'.format( weighted ) ]
	return outAr

def writeOutput( outFileStr, sampleMethylAr, sampleNamesAr, methTypes, info ):
	
	outFile = open( outFileStr, 'w' )
	# header
	outStr = info + '\nsample'
	for x in methTypes:
		outStr += '\tmC-'+x+'\ttC-'+x+'\twm-'+x
	outStr += '\n'
	outFile.write( outStr )
	
	# loop through samples
	for i in range(len(sampleNamesAr)):
		outStr = '{:s}\t{:s}\n'.format( sampleNamesAr[i], '\t'.join( sampleMethylAr[i] ) )
		outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	
	outPre = 'out'
	numProc = NUMPROC
	includeList = []
	excludeList = []
	sampleFile = False
	binTest = False
	isPrint = True
	methTypes = ['CG','CHG','CHH','C']
	startInd = 0
	
	for i in range( min(8, len(argv)-2) ):
		
		if argv[i] == '-f':
			sampleFile = True
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
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
		elif argv[i].startswith( '-i='):
			if len(excludeList) > 0:
				print( 'ERROR: cannot specify include and exclude lists' )
				exit()
			includeList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-x=' ):
			if len(includeList) > 0:
				print( 'ERROR: cannot specify include and exclude lists' )
				exit()
			excludeList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-m='):
			methTypes = argv[i][3:].split(',')
			startInd += 1
		elif argv[i] == '-b':
			binTest = True
			startInd += 1
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith('-'):
			print('ERROR: {:s} is not a valid parameter'.format( argv[i] ))
			exit()
	# end for

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
	processInputs( allCPath, sampleNamesAr, numProc, includeList, excludeList, methTypes, outPre, sampleFile, binTest, isPrint )

def printHelp():
	print ("Usage:\tpython3 genome_methylation.py [-q] [-f] [-b] [-p=num_proc]\n\t[-o=out_prefix] [-i=include_chrms | -x=exclude_chrms] [-m=meth_types] \n\t<allC_path> <sample_name> [sample_name]*")
	print('Compute genome wide for various samples')
	print('Optional:')
	print('-q\t\tquiet; do not print progress')
	print("-f\t\tsamples are listed in file (1 line per sample)")
	print("-b\t\tcompute weighted methylation for positions passing binomial test only")
	print("-i=include_chrms\tcomma-separated list of chmrs to ONLY include")
	print("\tor")
	print("-x=exclude_chrms\tcomma-separated list of chrms to exclude")
	print("\t\t[default includes all chrms in allC file]")
	print("-m=meth_types\tcomma-separated list methylation contexts to analyze [default CG,CHG,CHH,C]")
	print("-p=num_proc\tnumber of processors to use [default 1]")
	print("-o=out_prefix\tprefix to use for output file [default 'out']")

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
