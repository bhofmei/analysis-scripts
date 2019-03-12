import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python parse_pipeline_output.py [-s] [-o=out_id] [-u=underscore] <read_length> <fasta_index> <path_to_output>
# Folder with output files should contain both error and std out

def processInputs( pathRun, readLength, fastaFileStr, outID, underScores, simpleOutput ):
	# check that folder exists
	if os.path.exists( pathRun ) == False:
		print( 'ERROR: folder {:s} not found'.format( pathRun ) )
		exit()
	# initialize dictionary
	totalDict = {}
	# get run and error files
	runFiles = glob.glob( os.path.normpath( pathRun + '/*.o*'  ) )
	errorFiles = glob.glob( os.path.normpath( pathRun + '/*.e*'  ) )
	
	print( 'Reading', os.path.basename(fastaFileStr ) )
	cov = readIndex( fastaFileStr )	
	
	print( 'Reading output files' )
	for run in runFiles:
		n = getSampleName( run, underScores+1 )
		print( '-',n)
		if totalDict.get(n) == None:
			tempAr = readRunOutput( run )
			totalDict[n]=tempAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(n, run) )
	print( 'Reading error files' )
	for err in errorFiles:
		n = getSampleName( err, underScores+1 )
		print( '-',n)
		if totalDict.get(n) == None:
			print( 'ERROR: {:s} has error file but no run file'.format( n ) )
		else:
			tAr = totalDict.get(n)
			tmpAr = readErrorOutput( err )
			tmpAr2 = [tAr[0]] + tmpAr + tAr[1:] # input, forward once, forward multiple, ...
			totalDict[n] = tmpAr2
	
	print('Computing coverage')
	totalDict = computeCoverage(totalDict, cov, readLength)
	
	outFileStr = '{:s}methyl-stats.csv'.format( outID+'_' if outID != None else '' )
	if simpleOutput:
		outFileStr = outFileStr.replace('.csv', '_simple.csv' )
	print( 'Writing output to', outFileStr )
	writeOutput( totalDict, outFileStr, simpleOutput )
	print( 'Done.' )

def readIndex( fastaIndexStr ):
	baseCount = 0
	indexFile = open( fastaIndexStr, 'r' )
	for line in indexFile:
		lineAr = line.rstrip().split('\t')
		# (0) name (1) length (2) byte position (3) line length w/o \n
		# (4) line length w/ \n
		#if lineAr[0] not in ['mitochondria','chloroplast','ChrC','ChrM','lambda']:
		baseCount += int( lineAr[1] )
	indexFile.close()
	return baseCount

def getSampleName( fileStr, underScore ):
	baseName = os.path.basename( fileStr )
	fInd = 0
	for i in range( underScore ):
		tInd = baseName[fInd+1:].find( '_' )
		if tInd != -1:
			fInd += tInd + 1
	if fInd == 0:
		tInd = baseName.find( '.' )
		return baseName[:tInd]
	return baseName[:fInd]

def readRunOutput( runFileStr ):
	# (0) input reads (1) uniquely mapped (2) non-clonal
	# (3) percent remaining (4) non-conversion
	outAr = [-1]*5
	
	state = 0
	
	runFile = open( runFileStr, 'r' )
	
	for line in runFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		if state == 0 and lineAr[0] == 'There':
			# There are 199024702 total input reads
			outAr[0] = int( lineAr[2] )
			state = 1
		elif state == 1 and lineAr[0] == 'There':
			# There are 126515226 uniquely mapping reads, 63.5675997646 percent remaining
			outAr[1] = int( lineAr[2] )
			state = 2
		elif state == 2 and lineAr[0] == 'There':
			# There are 104890077 non-clonal reads, 52.7020394684 percent remaining
			outAr[2] = int( lineAr[2] )
			outAr[3] = float( lineAr[5] )
			state = 3
		elif state == 3 and lineAr[1] == 'non-conversion':
			# The non-conversion rate is 0.158129364499%
			outAr[4] = float(lineAr[4].replace("%",""))
			break
	# end for
	runFile.close()
	return outAr
		
def readErrorOutput( errorFileStr ):
	# (0) forward-mapped once rate, (1) forward-mapped multiple rate, (2) reverse-mapped once rate, (3) reverse-mapped multiple rate
	errorFile = open( errorFileStr, 'r' )
	outAr = [ -1 ] * 4
	# forward_one_alignment 
	state = 0
	
	for line in errorFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		# 199015159 (100.00%) were unpaired; of these:
		# 71720187 (36.04%) aligned 0 times
		# 68390844 (34.36%) aligned exactly 1 time
		# 58904128 (29.60%) aligned >1 times
		# 63.96% overall alignment rate
		if (state == 0 or state == 3) and lineAr[2] == 'aligned' :
			# aligned 0
			state += 1
		elif state == 1 and lineAr[2] == 'aligned':
			# forward aligned once
			outAr[0] = _cleanPercent(lineAr[1])
			state = 2
		elif state == 2 and lineAr[2] == 'aligned':
			# forward aligned multiple
			outAr[1] =  _cleanPercent(lineAr[1])
			state = 3
		elif state == 4 and lineAr[2] == 'aligned':
			# reverse aligned once
			outAr[2] = _cleanPercent(lineAr[1])
			state = 5
		elif state == 5 and lineAr[2] == 'aligned':
			# reverse aligned multiple
			outAr[3] =  _cleanPercent(lineAr[1])
			state = 6
		elif line.startswith( '[bam_sort_core]' ):
			break
	errorFile.close()
	return outAr

def _cleanPercent( inNumStr ):
	# inNum looks like (N%)
	# return N as float
	inNum = inNumStr.replace('(', '').replace(')', '').replace('%', '')
	try:
		return float(inNum)
	except ValueError:
		return -1

def computeCoverage( totalDict, genomeSize, readLength ):
	for samp in totalDict.keys():
		tmpAr = totalDict[samp]
		# (0) input reads (1) forward once (2) forward multiple (3) reverse once
		# (4) reverse multiple (5) unique mapped reads (6) non-clonal reads
		# (7) percent remaining (8) non-conversion
		# mappedBases = non-clonal reads * readLength
		mappedBases = tmpAr[6] * readLength
		cov = mappedBases / genomeSize
		totalDict[samp] = [readLength] + tmpAr + [cov]
	# end for
	return totalDict

def writeOutput( totalDict, outFileStr, simpleOutput ):
	outFile = open( outFileStr, 'w' )
	if simpleOutput:
		headerStr = 'no,sample,read length,input reads,non-conversion(%),mapped reads(unique non-clonal),mapped reads(%),genome coverage' 
	else:
		headerStr = 'sample,input reads,forward mapped once(%),forward mapped multiple(%),reverse mapped once(%),reverse mapped multiple(%),uniquely mapped reads,non-clonal reads,reads remaining(%),non-conversion(%),genome coverage'
	outFile.write( headerStr + '\n' )
	count = 1
	for sample in sorted(totalDict.keys()):
		tup = totalDict.get( sample )
		# (0) read length (1) input reads (2) forward once (3) forward multiple 
		# (4) reverse once (5) reverse multiple (6) unique mapped reads 
		# (7) non-clonal reads (8) percent remaining (9) non-conversion 
		# (10) coverage
		ar = [ '{:d}'.format(x) for x in tup[0:2] ]
		ar += [ '{:.2f}'.format(x) for x in tup[2:6] ]
		ar += [ '{:d}'.format(x) for x in tup[6:8] ]
		ar += [ '{:.2f}'.format(tup[8]) ]
		ar += [ '{:.3f}'.format(tup[9]) ]
		ar += [ '{:.2f}'.format(tup[10]) ]
		
		if simpleOutput:
			newAr = [ar[0]] + [ar[1]] + [ar[9]] + [ar[7]] + [ar[8]] + [ar[10]]
			outStr = str(count) + ',' + sample + ',' + ','.join(newAr) + '\n'
			count += 1
		else:
			outStr = sample+','+','.join(ar[1:])+'\n'

		outFile.write( outStr )
	outFile.close()
	
def parseInputs( argv ):
	outID = None
	underScores = 0
	simpleOutput = False
	startInd = 0
	
	for i in range(min(3,len(argv)-3)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-s':
			simpleOutput = True
			startInd += 1
		elif argv[i].startswith( '-u=' ):
			try:
				underScores = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of underscores must be an integer' )
				exit()
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	try:
		readLength = int( argv[startInd] )
	except ValueError:
		print( 'ERROR: read length must be integer' )
		exit()
	
	fastaFileStr = argv[startInd+1]
	pathRun = argv[startInd+2]

	processInputs( pathRun, readLength, fastaFileStr, outID, underScores, simpleOutput )

def printHelp():
	print ("Usage: python parse_pipeline_output.py [-h] [-s] [-o=out_id] [-u=underscore] <read_length> <fasta_index> <path_to_output>")
	print( 'Required: ')
	print( 'read_length\tapprox read length before trimming')
	print( 'fasta_index\tgenome fasta index (chrm name and size)')
	print( 'path_to_run_output\tpath to *.o* and *.e*' )
	print( 'Optional: ' )
	print( '-s\t\tsimpler output format; requires -f' )
	print( '-u=underscore_count\tunderscores in sample id name [default 0]' )
	print( '-o=out_id\toutput identifier' )
	print( '-h\t\tprint this help and exit' )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
