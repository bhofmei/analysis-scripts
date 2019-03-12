import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python parse_smrna_output.py [-s] [-o=out_id] [-u=underscore] <path_to_output>
# Folder with output files should contain both error and std out

def processInputs( pathRun, outID, underScores, simpleOutput ):
	# check that folder exists
	if os.path.exists( pathRun ) == False:
		print( 'ERROR: folder {:s} not found'.format( pathRun ) )
		exit()
	# initialize dictionary
	totalDict = {}
	# get run and error files
	runFiles = glob.glob( os.path.normpath( pathRun + '/*.o*'  ) )
	errorFiles = glob.glob( os.path.normpath( pathRun + '/*.e*'  ) )
	
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
			#tmpAr2 = tAr[0:4] + tmpAr + tAr[4:]
			tmpAr2 = tAr + tmpAr
			totalDict[n] = tmpAr2
	
	outFileStr = 'smrna_output{:s}.csv'.format( '_'+outID if outID != None else '' )
	if simpleOutput:
		outFileStr = outFileStr.replace('.csv', '_simple.csv' )
	print( 'Writing output to', outFileStr )
	writeOutput( totalDict, outFileStr, simpleOutput )
	print( 'Done.' )

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
	# (0) input reads (1) input bases (2) reads with adapter (3) surviving bases
	# (4) 3` adapter trimmed (5) 5` adapter trimmed (6) surviving reads
	outAr = [0]*7
	
	state = 0
	
	runFile = open( runFileStr, 'r' )
	
	for line in runFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		if state == 0:
			#print( lineAr )
			if lineAr[0] == 'Total' and lineAr[1] == 'reads':
				# Total reads processed: 7,848,544
				outAr[0] += int( lineAr[3].replace( ',', '' ) )
			elif lineAr[0] == 'Reads' and lineAr[2] == 'adapters:':
				# Reads with adapters: 7,531,380 (96.0%)
				outAr[2] += int( lineAr[3].replace( ',', '' ) )
			elif lineAr[0] == 'Reads' and lineAr[1] == 'written':
				# Reads written (passing filters):    12,039,362 (100.0%)
				outAr[6] += int( lineAr[4].replace( ',', '' ) )
			elif lineAr[0] == 'Total' and lineAr[1] == 'basepairs':
				# Total basepairs processed: 1,811,125,242 bp
				outAr[1] += int( lineAr[3].replace(',', '' ) )
			elif lineAr[0] == 'Total' and lineAr[1] == 'written':
				# Total written (filtered):  1,810,609,728 bp (100.0%)
				outAr[3] += int( lineAr[3].replace( ',', '' ) )
			elif lineAr[0] == "===" and lineAr[1] == 'Adapter':
				state = 1
		elif state == 1 and lineAr[0] == 'Sequence:':
			# Sequence: TGGAATTCTCGGGTGCCAAGGAACTCCAGT; Type: regular 3'; Length: 30; Trimmed: 7524980 times.
			
			outAr[4] = int( lineAr[-2] )
			state = 2
		elif state == 2 and lineAr[0] == 'Sequence:':
			# Sequence: TGTGTTCTCAGGTCGCCCCTG; Type: regular 5'; Length: 21; Trimmed: 6400 times.
			outAr[5] = int( lineAr[-2] )
			break
	# end for
	runFile.close()
	
	# percent remaining
	outAr += [ float( outAr[6] ) / float(outAr[0] ) * 100 ]
	# approx read length
	approxReadLength = math.floor( float(outAr[1]) / float(outAr[0]) )
	outAr2 = [approxReadLength] + outAr
	print( outAr2 )
	return outAr2
		
def readErrorOutput( errorFileStr ):
	# (0) mapped reads
	errorFile = open( errorFileStr, 'r' )
	outAr = [ -1 ] * 2
	# (0) reads mapped (1) percent mapped 
	state = 0
	
	for line in errorFile:
		line = line.rstrip()
		if line.startswith( '# reads with at least one reported alignment' ):
			## reads with at least one reported alignment: 593022 (14.35%)
			lineAr = line.split(' ')
			p = lineAr[-1].replace('(','').replace(')','').replace('%','')
			outAr[0] = int( lineAr[-2] )
			outAr[1] = float( p )
			#state += 1
		elif line.startswith( '[bam_sort_core]' ):
			break
	errorFile.close()
	print( outAr )
	return outAr

def writeOutput( totalDict, outFileStr, simpleOutput ):
	# totalDict[n]=run_array
	outFile = open( outFileStr, 'w' )
	if simpleOutput:
		headerStr = 'no,sample,read length,input reads,reads after adapter trimming,percent remaining,mapped reads,percent mapped' 
	else:
		headerStr = 'sample,read length,input reads,input bases,reads with adapters,surviving bases,3 adapter trimmed,5 adapter trimmed,surviving reads,percent remaining,reads mapped,percent mapped'
	outFile.write( headerStr + '\n' )
	#print( totalDict['line1-G3-rep1'])
	#print( sorted(totalDict.keys()) )
	count = 1
	for sample in sorted(totalDict.keys()):
		tup = totalDict.get( sample )
		# (0) read length (1) input reads (2) input bases (3) reads with adapter 
		# (4) surviving bases (5) 3` adapter trimmed (6) 5` adapter trimmed 
		# (7) surviving reads(8) percent remaining (9) mapped reads 
		# (10) percent mapped
		ar = [ '{:d}'.format(x) for x in tup[0:8] ]
		ar += [ '{:.2f}'.format(tup[8] ) ]
		ar += [ '{:d}'.format(tup[9])]
		ar += [ '{:.2f}'.format(tup[10]) ]
		
		if simpleOutput:
			newAr = [ar[0]] + [ar[1]] + [ar[7]] + [ar[8]] + [ar[9]] + [ar[10]]
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
	underScores = 0
	simpleOutput = False
	startInd = 0
	
	for i in range(min(3,len(argv)-1)):
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
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for 
	pathRun = argv[startInd]

	processInputs( pathRun, outID, underScores, simpleOutput )

def printHelp():
	print ("Usage: python parse_smrna_output.py [-s] [-o=out_id] [-u=underscore] <path_to_output>")
	print( 'Required: ')
	print( 'path_to_run_output\tpath to *.o* and *.e*' )
	print( 'Optional: ' )
	print( '-s\t\tsimpler output format' )
	print( '-u=underscore_count\tunderscores in sample id name [default 0]' )
	print( '-o=out_id\toutput identifier' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
