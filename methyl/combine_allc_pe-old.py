import sys, multiprocessing, subprocess, os, gzip

# Usage: combine_allc_pe.py [-f] [-p=num_proc] [-o=out_id] [-c=chrm_list | -cf=fasta_index] <allc_path> <sample_name> [sample_name]*
# combine many allc files at basepair level into one allc file

NUMPROC=1
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def processInputs( allCPath, sampleNamesAr, numProc, chrmList, fastaIndex, outID, sampleFile, isPrint ):

	
	if outID == None:
		outID = 'combined'
	
	if sampleFile:
		if isPrint:
			print( 'Reading file of samples names...' )
		sampleNamesAr = readSampleFile( sampleNamesAr[0] )
	if fastaIndex != None:
		chrmList = readFastaIndex( fastaIndex )
	
	info = '#from_script: combine_allc_pe.py; samples: {:s}'.format( ','.join(sampleNamesAr) )

	#if isPrint:
		#print( '\nChromosomes: {:s}\nSamples included: {:s}\n'.format( ' '.join(chrmList), ' '.join(sampleNamesAr) ) )
	
	# check for all allC files before doing any work
	#for sample in sampleNamesAr:
	#	if checkFiles( allCPath, sample, chrmList ) == False:
	#		exit()
	
	#  loop through chromosomes
	if isPrint:
		print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processChrm, args=(sampleNamesAr, allCPath, chrm, outID, info, isPrint) ) for chrm in chrmList ]
	suc = [ p.get() for p in results ]
	if isPrint:
		print( 'Done' )

def readSampleFile( fileStr ):
	
	if os.path.isfile(fileStr) == False:
		print( 'ERROR:', os.path.basename(fileStr), 'does not exist' )
		exit()
	sampleAr = []
	inFile = open( fileStr, 'r' )
	for line in inFile:
		name = line.rstrip()
		sampleAr += [ name ]
	inFile.close()
	return sampleAr

def readFastaIndex( fastaIndexStr ):
	if os.path.isfile(fastaIndexStr) == False:
		print( 'ERROR:', os.path.basename(fastaIndexStr), 'does not exist' )
		exit()
	chrmList = []
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		chrmList += [ lineAr[0] ]
	fastaIndex.close()
	return chrmList

def checkFile( allcPath, sample, chrm ):
	
	if os.path.isfile( os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sample, chrm ) ) ) == False and os.path.isfile( os.path.normpath('{:s}/allc_{:s}_{:s}.tsv.gz'.format( allcPath, sample, chrm ) ) ) == False:
		return False
	return True

def processChrm( sampleNamesAr, allcPath, chrm, outID, info, isPrint ):
	# initialize dictionary
	methDict = {}
	
	# loop through samples
	for sample in sampleNamesAr:
		if checkFile(allcPath, sample, chrm):
			allCFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sample, chrm ) )
			if isPrint:
				print( 'Reading {:s}'.format( allCFileStr ) )
			methDict = readAllC( allCFileStr, methDict )
	
	# write output
	if len(methDict.keys()) > 0:
		outFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, outID, chrm ) )
		if isPrint:
			print( 'Writing', outFileStr )
		writeOutput( outFileStr, methDict, chrm, info )
	
	# delete dict
	del( methDict )
	return True
	
def readAllC( allCFileStr, inDict ):
	allCFile = openFile( allCFileStr )
	# headers
	for line in allCFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		pos = int( lineAr[1] )
		
		tup = inDict.get( pos )
		# doesn't exist in dict -> add
		if tup == None:
			inDict[pos] = [ lineAr[2], lineAr[3], int( lineAr[4] ), int( lineAr[5] ) ]
		else:
			# check strand and type match
			if lineAr[2] != tup[0]:
				print( 'WARNING: strand error for {:s} at position {:d}'.format( os.path.basename( allCFileStr ), pos ) )
			elif lineAr[3] != tup[1]:
				# count N's
				cn1 = tup[1].count('N')
				cn2 = lineAr[3].count('N')
				if cn1 == 0 and cn2 == 0:
					print( 'WARNING: methylation type error for {:s} at position {:d}'.format( os.path.basename( allCFileStr ), pos ) )
				# replace type if new file has fewer N's
				if cn2 < cn1:
					inDict[pos][1] = lineAr[3]
			inDict[pos][2] += int( lineAr[4] )
			inDict[pos][3] += int( lineAr[5] )
		# end else
	# end for
	allCFile.close()
	return inDict

def openFile( allcFileStr ):
	if os.path.isfile( allcFileStr ) == False and os.path.isfile( allcFileStr + '.gz' ):
		return gzip.open( allcFileStr+'.gz', 'rt' )
	else:
		return open( allcFileStr, 'r' )
			
def writeOutput( outFileStr, methDict, chrm, info ):
	# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
	# (6) methylated
	
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' )
	
	# loop through positions
	for pos in sorted( methDict.keys() ):
		info = methDict[pos]
		outFile.write( '{:s}\t{:d}\t{:s}\t{:s}\t{:d}\t{:d}\t1\n'.format( chrm, pos, info[0], info[1], info[2], info[3] ) )
	# end for
	outFile.close()

def parseInputs( argv ):
	outID = None
	numProc = NUMPROC
	sampleFile = False
	chrmList = None
	fastaIndex = None
	isPrint = True
	startInd = 0
	
	for i in range(min(5,len(argv)-2)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-f':
			sampleFile = True
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
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )

			except ValueError:
				print( 'WARNING: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
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
	
	processInputs( allCPath, sampleNamesAr, numProc, chrmList, fastaIndex, outID, sampleFile, isPrint )

def printHelp():
	
	print( 'Usage:\tpython combine_allc_pe.py [-q] [-h] [-f] [-p=num_proc] [-o=out_id]\n\t[-c=chrm_list | -cf=fasta_index] <allc_path> <sample_name> [sample_name]*' )
	print()
	print( 'Required:' )
	print( 'allc_path\tpath to allC files' )
	print( 'sample_name\tname of sample; used to find allC files\n\t\twhen "-f" flag set, file with sample names listed one per line' )
	print()
	print( 'Optional:' )
	print( '-q\t\tquiet; do not print progress' )
	print( '-h\t\tprint help message and exit' )
	print( '-f\t\tsample names are in file' )
	print( '-p=num_proc\tnumber of processors to use [default {:d}]'.format( NUMPROC) )
	print( '-o=out_id\toutput file identifier [default "combined"]' )
	print( '-c=chrm_list\tcomma-separated list of chrms to use [default for arabidopsis]' )
	print( '-cf=fasta_index\tfasta index file with chrms to use' )
	
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
