import sys, math, glob, multiprocessing, subprocess, os, pickle, gzip

# Usage: python tabulate_methyl_allc_pe.py [-q] [-h] [-z] [-v=min_cov] [-m=meth_type] [-o=out_id] [-p=num_proc] <allc_path> <sample1> [sampleN]*

MINCOV=1
NUMPROC=1
ZFILL=12

def processInputs( allcPath, samplesAr, outId, minCov, numProc, methType, isCompress, isPrint ):
	
	if isPrint:
		print('Min cov:', minCov)
		print('Methyl type:', methType)
		print()
		print('Begin processing', len(samplesAr), 'samples with', numProc, 'processors')
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(processSample, args=(allcPath, sample, minCov, methType, isPrint) ) for sample in samplesAr ]
	chrmList = [ p.get() for p in results ]
	
	# update sample list
	updateSamples = []
	for i in range(len(samplesAr)):
		updateSamples += [samplesAr[i]] if chrmList[i] != False else []
	
	allChrmList = mergeChrms( chrmList )
	
	if isPrint:
		print('Begin merging', len(allChrmList), 'chromosomes for', len(updateSamples), 'samples with', numProc, 'processors' )
	
	results = [ pool.apply_async(processChrm, args=(chrm, updateSamples, isPrint) ) for chrm in allChrmList ]
	savedChrms = [ p.get() for p in results ]
	
	# final merge
	savedChrms.sort()
	
	outFileStr = '{:s}_tabulated_{:s}.tsv{:s}'.format(outId, methType.lower(), ('.gz' if isCompress else '') )
	info = '#from_script: tabulate_methyl_allc_pe.py; min_cov: {:d}; methyl_type: {:s}; samples: {:s}'.format( minCov, methType, ','.join(updateSamples) )
	
	if isPrint:
		print( 'Writing output to', outFileStr )
	
	processOutput( outFileStr, savedChrms, updateSamples, info, isCompress )
	
	if isPrint:
		print( 'Done' )

def processSample( allcPath, sampleName, minCov, methType, isPrint ):
	if isPrint:
		print('~Processing', sampleName)
	
	allcFileStr = getAllcFileName( allcPath, sampleName )
	
	if allcFileStr == False:
		return False
	else:
		return readAllc( allcFileStr, minCov, methType, sampleName, isPrint )

def getAllcFileName(allcPath, sampleName):
	allcFileStr = '{:s}/allc_{:s}.tsv'.format(allcPath, sampleName)
	
	if os.path.isfile(allcFileStr):
		return allcFileStr
	elif os.path.isfile( allcFileStr + '.gz' ):
		return allcFileStr + '.gz'
	else:
		print('ERROR: No allc file for {:s}...skipping'.format(sampleName) )
		return False
	
def readAllc( allcFileStr, minCov, methType, sampleName, isPrint ):
	
	tmpFilePre = 'tmp_'+sampleName
	allcFile = openAllc( allcFileStr )
	
	currChrm = None
	outDict = {}
	chrmAr = []
	
	for line in allcFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if line.startswith('#') or len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		chrm = lineAr[0]
		pos = lineAr[1]
		strand = lineAr[2]
		cov = int(lineAr[5])
		mType = decodeMethType ( lineAr[3] )
		isMeth = lineAr[6]
		
		# new chrm - save results
		if chrm != currChrm:
			if isPrint and currChrm != None:
				print( '~~Write', currChrm )
			r = savePickle( tmpFilePre, currChrm, outDict )
			chrmAr += [r] if r != '' else []
			outDict = {}
			currChrm = chrm
		
		# do we save this
		if (cov >= minCov) and (methType == 'C' or mType == methType):
			outName = pos.zfill(ZFILL) + '.' + strand
			outDict[outName] = isMeth
	# end for line
	
	# save the last dict
	if isPrint:
		print( '~~Write', currChrm )
	r = savePickle( tmpFilePre, currChrm, outDict )
	chrmAr += [r] if r != '' else []
	
	allcFile.close()
	return chrmAr
	
def openAllc( allcFileStr ):
	if allcFileStr.endswith('.gz'):
		return gzip.open( allcFileStr, 'rt' )
	else:
		return open( allcFileStr, 'r' )

def decodeMethType( mStr ):

	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return 'CNN'
	else:
		return 'CHH'

def mergeChrms( chrmLists ):
	
	chrmSet = set()
	# loop thru samples
	for sList in chrmLists:
		# loop thru chrms
		for chrm in sList:
			if chrm != '':
				chrmSet.add( chrm )
	
	return sorted(list(chrmSet))
		
def processChrm( chrm, samplesAr, isPrint ):
	if isPrint:
		print('~Merging', chrm)
		
	tmpFilePre = 'tmp'
	nSamples = len(samplesAr)
	outDict = {}
	
	# loop through samples
	for i in range(nSamples):
		sampleName = samplesAr[i]
		if isPrint:
			print('~~Read', sampleName)
		sampleFilePre = 'tmp_' + sampleName
		tmpDict = openPickle( sampleFilePre, chrm )
		# if we have tmpDict, merge and remove pickle file
		if tmpDict != None:
			outDict = mergeDict( outDict, tmpDict, nSamples, i)
			removePickle( sampleFilePre, chrm )
	# end for i
	
	# save to pickle
	chrmSuccess = savePickle( tmpFilePre, chrm, outDict )
	return chrmSuccess
			
def mergeDict(origDict, newDict, n, i):
	# loop thru newDict
	for key in newDict.keys():
		if origDict.get(key) == None:
			# not in dict, must add
			origDict[key] = ['NA']*n
		origDict[key][i] = newDict[key] # save
	# end for key
	return origDict
	
def processOutput( outFileStr, savedChrms, updateSamples, info, isCompress ):
	
	outFile = openOutFile( outFileStr, isCompress )
	headerAr = ['chrm', 'pos', 'strand'] + updateSamples
	
	outFile.write(info+'\n')
	outFile.write('\t'.join(headerAr) + '\n')
	
	tmpFilePre = 'tmp'
	
	# loop thru chrms
	for chrm in savedChrms:
		tmpDict = openPickle( tmpFilePre, chrm )
		if tmpDict != None:
			outStr = formatOutputDict( chrm, tmpDict )
			outFile.write(outStr)
			removePickle( tmpFilePre, chrm )
	# end for chrm
	outFile.close()
	return True
	
def formatOutputDict( chrm, inDict ):
	outStr = ''
	for outName in sorted(inDict.keys()):
		vals = inDict[outName]
		pos, strand = outName.split('.')
		outStr += '{:s}\t{:d}\t{:s}\t'.format( chrm, int(pos), strand ) + '\t'.join(vals) + '\n'
	return outStr
	
def openOutFile( outFileStr, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	return outFile

def savePickle( filePrefix, chrm, outDict ):
	if len(outDict.keys()) == 0:
		# empty dict - return
		return ''
	else:
		# pickle result
		outFileStr = filePrefix + '_' + chrm + '.pick'
		pickle.dump( outDict, open( outFileStr, 'wb'), protocol=3 )
		return chrm

def openPickle( filePrefix, chrm ):
	outFileStr = filePrefix + '_' + chrm + '.pick'
	if os.path.exists( outFileStr ):
		outDict = pickle.load( open(outFileStr, 'rb') )
		return outDict
	else:
		return None

def removePickle( filePrefix, chrm ):
	outFileStr = filePrefix + '_' + chrm + '.pick'
	os.remove( outFileStr )

def parseInputs( argv ):

	outId = 'out'
	minCov = MINCOV
	numProc = NUMPROC
	methType = 'C'
	isPrint = True
	isCompress = False
	startInd = 0
	
	for i in range(min(6, len(argv))):
		if argv[i] =='-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i].startswith('-v='):
			try:
				mCov =int(argv[i][3:])
				if (mCov < 1):
					print('WARNING: Min cov must be at least 1...using default', MINCOV)
				else:
					minCov = mCov
			except ValueError:
				print('WARNING: Min cov must be integer...using default', MINCOV)
			startInd += 1
		elif argv[i].startswith('-p='):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'WARNING: Num processors must be integer...using default', NUMPROC )
			startInd += 1
		elif argv[i].startswith('-m='):
			mType = argv[i][3:].upper()
			if mType not in ['C', 'CG', 'CHG', 'CHH', 'CNN']:
				print('ERROR: Invalid methylation type. Choose one of "C", "CG", "CHG", "CHH"')
				exit()
			else:
				methType = 'C' if mType == 'CNN' else mType
				startInd += 1
		elif argv[i].startswith('-o='):
			outId = argv[i][3:]
			startInd += 1
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith('-'):
			print('ERROR: Invalid parameter', argv[i])
			exit()
	# end for i
	
	allcPath = argv[startInd]
	samplesAr = argv[startInd+1:]
			
	processInputs( allcPath, samplesAr, outId, minCov, numProc, methType, isCompress, isPrint )

def printHelp():
	print( 'Usage: python\ttabulate_methyl_allc_pe.py [-q] [-h] [-z]' )
	print( '\t\t[-v=min_cov] [-m=meth_type] [-o=out_id] [-p=num_proc]' )
	print( '\t\t<allc_path> <sample1> [sampleN]*' )
	print()
	print( 'Required:' )
	print( 'allc_path\tdirectory with allC files; 1 file per sample')
	print( 'sampleN\t\tsamples to include' )
	print()
	print( 'Optional:' )
	print( '-h\t\tprint help message and exit' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-z\t\tcompress output file with gzip' )
	print( '-v=min_cov\tmin # reads per position for it to be included\n\t\t[default {:d}, i.e. include all]' )
	print( '-m=meth_type\tmethylation type to include; must be one of: "CG",\n\t\t"CHG", "CHH", or "C" [default C]' )
	print( '-o=out_id\toutput file prefix [default "out"]')
	print( '-p=num_proc\tnumber of processors to use [default {:d}]'.format(NUMPROC) )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
