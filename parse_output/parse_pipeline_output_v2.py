import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python parse_pipeline_output.py [-s] [-f=fasta_file | fasta_index] [-o=out_id] [-u=underscore] <path_to_output>
# Folder with output files should contain both error and std out

def processInputs( pathRun, fastaFileStr, outID, underScores, simpleOutput ):
	# check that folder exists
	if os.path.exists( pathRun ) == False:
		print( 'ERROR: folder {:s} not found'.format( pathRun ) )
		exit()
	# initialize dictionary
	totalDict = {}
	# get run and error files
	runFiles = glob.glob( os.path.normpath( pathRun + '/*.o*'  ) )
	errorFiles = glob.glob( os.path.normpath( pathRun + '/*.e*'  ) )
	
	cov = -1
	# get fasta or fasta index if specified
	if fastaFileStr != None:
		print( 'Reading', os.path.basename(fastaFileStr ) )
		if fastaFileStr.endswith( '.fai' ):
			cov = readIndex( fastaFileStr )	
		elif fastaFileStr != None:
			cov = readFasta( fastaFileStr )
		print( 'Genome size: {:.4f} MB'.format( cov / 1000000.0 ) )
	
	print( 'Reading output files for...' )
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
			tmpAr2 = tAr[0:4] + tmpAr + tAr[4:]
			totalDict[n] = tmpAr2
	
	outFileStr = 'methyl_output{:s}.csv'.format( '_'+outID if outID != None else '' )
	if simpleOutput:
		outFileStr = outFileStr.replace('.csv', '_simple.csv' )
	print( 'Writing output to', outFileStr )
	writeOutput( totalDict, cov, outFileStr, simpleOutput )
	print( 'Done.' )

def readIndex( fastaIndexStr ):
	baseCount = 0
	indexFile = open( fastaIndexStr, 'r' )
	for line in indexFile:
		lineAr = line.rstrip().split('\t')
		# (0) name (1) length (2) byte position (3) line length w/o \n
		# (4) line length w/ \n
		if lineAr[0] not in ['mitochondria','chloroplast','ChrC','ChrM','lambda']:
			baseCount += int( lineAr[1] )
	indexFile.close()
	return baseCount

def readFasta( fastaFileStr ):
	baseCount = 0
	fastaFile = open( fastaFileStr, 'r' )
	include = False
	for line in fastaFile:
		line = line.rstrip()
		# chromosome headers
		if line.startswith('>'):
			# check chloroplast/mitochondria
			if line.startswith( '>chloroplast' ) or line.startswith( '>mitochondria'):
				include = False
			else:
				include = True
		# sequence to be included
		elif include:
			baseCount += len( line )
	fastaFile.close()
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
	# (0) input reads (1) input bases (2) surviving reads  
	# (3) quality trimmed bases (4) uniquely mapped (5) non-clonal
	# (6) percent remaining (7) genome-pvalue (8) non-conversion
	outAr = [0]*4 + [-1] * 5
	
	state = 0
	
	runFile = open( runFileStr, 'r' )
	
	for line in runFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		if state == 0 and lineAr[0] == 'There':
			# There are 12039545 total input reads
			outAr[0] = int( lineAr[2] )
			state = 1
		elif state == 1:
			if lineAr[0] == 'Reads' and lineAr[1] == 'written':
				# Reads written (passing filters):    12,039,362 (100.0%)
				outAr[2] += int( lineAr[4].replace( ',', '' ) )
			elif lineAr[0] == 'Total' and lineAr[1] == 'basepairs':
				# Total basepairs processed: 1,811,125,242 bp
				outAr[1] += int( lineAr[3].replace(',', '' ) )
			elif lineAr[0] == 'Total' and lineAr[1] == 'written':
				# Total written (filtered):  1,810,609,728 bp (100.0%)
				outAr[3] += int( lineAr[3].replace( ',', '' ) )
			elif lineAr[0] == "Begin" and lineAr[1] == "converting":
				state = 2
		elif state == 2 and lineAr[0] == 'There':
			# There are 5202567 uniquely mapping reads, 43.2123223926 percent remaining
			outAr[4] = int( lineAr[2] )
			state = 3
		elif state == 3 and lineAr[0] == 'There':
			# There are 4678906 non-clonal reads, 38.8628141678 percent remaining
			outAr[5] = int( lineAr[2] )
			outAr[6] = float( lineAr[5] )
			state = 4
		elif state == 4 and lineAr[2] == "p-value":
			# The minimum p-value in the control genome is 0.00159876921865
			outAr[7] = float( lineAr[8] )
			state = 5
		elif state == 5:
			# The non-conversion rate is 0.293098889696%
			outAr[8] = float(lineAr[4].replace("%",""))
			break
	# end for
	runFile.close()
	return outAr
		
def readErrorOutput( errorFileStr ):
	# (0) forward-mapped rate (1) reverse-mapped rate
	errorFile = open( errorFileStr, 'r' )
	outAr = [ -1 ] * 2
	# forward_one_alignment 
	state = 0
	
	for line in errorFile:
		line = line.rstrip()
		if line.startswith( '# reads with at least one reported alignment' ):
			lineAr = line.split(' ')
			p = lineAr[-1].replace('(','').replace(')','').replace('%','')
			outAr[state] = float( p )
			state += 1
		elif line.startswith( '[bam_sort_core]' ):
			break
	errorFile.close()
	return outAr

def writeOutput( totalDict, cov, outFileStr, simpleOutput ):
	# totalDict[n]=run_array
	outFile = open( outFileStr, 'w' )
	if simpleOutput:
		headerStr = 'no,sample,read length,input reads,non-conversion%,mapped reads,percent mapped reads,genome coverage' 
	else:
		headerStr = 'sample,input reads,input bases,surving reads,quality trimmed bases,percent forward mapped,percent reverse mapped,uniquely mapped reads,non-clonal reads,percent reads remaining, genome pvalue,non-conversion percent'
		if cov != -1:
			headerStr += ',coverage,read length'
	outFile.write( headerStr + '\n' )
	#print( totalDict['line1-G3-rep1'])
	#print( sorted(totalDict.keys()) )
	count = 1
	for sample in sorted(totalDict.keys()):
		tup = totalDict.get( sample )
		# (0) input reads (1) input bases (2) surviving reads  
		# (3) quality trimmed bases (4) forward-mapped rate 
		# (5) reverse-mapped rate (6) uniquely mapped (7) non-clonal
		# (8) percent remaining (9) genome-pvalue (10) non-conversion
		ar = [ '{:d}'.format(x) for x in tup[0:4] ]
		ar += [ '{:.2f}'.format(x) for x in tup[4:6] ]
		ar += [ '{:d}'.format(x) for x in tup[6:8]]
		ar += [ '{:.5f}'.format(x) for x in tup[8:] ]
		if cov != -1:
			# mappedBases = uniquely mapped reads * floor( processBases / processReads )
			approxReadLength = math.floor( float(tup[1]) / float(tup[0]) )
			mappedBases = tup[6] * approxReadLength
			ar += [ '{:.3f}'.format( mappedBases/cov) ] + [ '{:d}'.format( int( approxReadLength) ) ]
		
		if simpleOutput:
			newAr = [ar[-1]] + [ar[0]] + [ar[10]] + [ar[7]] + [ar[8]] + [ar[-2]]
			outStr = str(count) + ',' + sample + ',' + ','.join(newAr) + '\n'
			count += 1
		else:
			outStr = sample+','+','.join(ar)+'\n'
		#except TypeError:
		#	print( sample )
			#print( tup )
			#print( ar )
			#print( fdrAr )
		outFile.write( outStr )
	outFile.close()
	
def parseInputs( argv ):
	outID = None
	fastaFileStr = None
	underScores = 0
	simpleOutput = False
	startInd = 0
	
	for i in range(min(3,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			fastaFileStr = argv[i][3:]
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
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for 
	if simpleOutput and fastaFileStr == None:
		print( 'ERROR: simple output format requires fasta file' )
		exit()
	pathRun = argv[startInd]

	processInputs( pathRun, fastaFileStr, outID, underScores, simpleOutput )

def printHelp():
	print ("Usage: python parse_pipeline_output.py [-s] [-f=fasta_file | fasta_index] [-o=out_id] [-u=underscore] <path_to_output>")
	print( 'Required: ')
	print( 'path_to_run_output\tpath to *.o* and *.e*' )
	print( 'Optional: ' )
	print( '-s\t\tsimpler output format; requires -f' )
	print( '-f=fasta\tfasta genome file or fasta index[default none]\n\t\tused to compute coverage' )
	print( '-u=underscore_count\tunderscores in sample id name [default 0]' )
	print( '-o=out_id\toutput identifier' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
