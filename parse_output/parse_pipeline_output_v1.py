import sys, os, glob, math

# Usage: python parse_pipeline_output.py [-f=fasta_file | fasta_index] [-o=out_id] [-u=underscore] <path_to_run_output> <path_to_correction_output>
# parse run/script information from output/error files

UNDERSCORE=1

def processInputs( pathRun, pathCorrection, fastaFileStr, outID, underScores ):
	uPathRun = checkFolder( pathRun )
	uPathCorrection = checkFolder( pathCorrection )
	if uPathRun == False or uPathCorrection == False:
		print( 'Error: folder {:s} or {:s} not found'.format( pathRun, pathCorrection ) )
	totalDict = {}
	runFiles = glob.glob(uPathRun+'*.py.o*')
	errorFiles = glob.glob(uPathRun+'*.py.e*')
	corFiles = glob.glob(uPathCorrection+'*.sh.o*')
	
	cov = -1
	if fastaFileStr != None:
		print( 'Reading', fastaFileStr )
		if fastaFileStr.endswith( '.fai' ):
			cov = readIndex( fastaFileStr )	
		elif fastaFileStr != None:
			cov = readFasta( fastaFileStr )
	print( 'Reading run files for...')
	for run in runFiles:
		n = getSampleName( run, underScores+1 )
		print( '-',n)
		if totalDict.get(n) == None:
			tempAr = readRunOutput( run )
			totalDict[n]=tempAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(n, run) )
	print( 'Reading error files for...' )
	for err in errorFiles:
		n = getSampleName( err, underScores+1 )
		print( '-',n)
		if totalDict.get(n) == None:
			print( 'ERROR: {:s} has error file but no run file'.format( n ) )
		else:
			tAr = totalDict.get(n)
			tmpAr = readRunError( err )
			tmpAr2 = tAr[0:4] + tmpAr + tAr[4:]
			totalDict[n] = tmpAr2
	print( 'Reading correction files for...' )
	for cor in corFiles:
		n = getSampleName( cor, underScores+1 )
		print('-',n )
		if totalDict.get(n) == None:
			print( 'ERROR: {:s} has correction file but not run file'.format(n ))
		else:
			tAr = totalDict.get(n)
			nonConversion, tempDict = readCorrentionOutput( cor )
			if tAr[-1] == -1:
				tAr[-1] = nonConversion
			#print( tAr )
			totalDict[n] = ( tAr, tempDict )
	#print( totalDict )
	# write output
	print( 'Writing output...')
	outFileStr = 'run_correction_output{:s}.csv'.format( '_'+outID if outID != None else '' )
	writeOutput( totalDict, cov, outFileStr )
	print( 'Done.' )

def checkFolder( inputStr ):
	if os.path.exists( inputStr ) == False:
		return False
	if inputStr.endswith( '/' ):
		return inputStr
	return inputStr + '/'
	
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
	baseName = baseName.replace( 'run_info_', '' )
	fInd = 0
	for i in range( underScore ):
		tInd = baseName[fInd+1:].find( '_' )
		if tInd != -1:
			fInd += tInd + 1
	if fInd ==0:
		tInd = baseName.find( '.' )
		return baseName[:tInd]
	return baseName[:fInd]

		
def readRunOutput( runFileStr ):
	# input_reads trimmed_reads x4 adapter_trimmed x4 quality_trimmed x4
	# uniquely_mapped non_clonal % remaining non_conversion_rate
	
	# [processed_reads, processed_bases, adapter_trimmed_reads, quality_trimmed_bases, unique, non_clonal, remaining, pvalue, nonConversion]
	outAr = [0]*4 + [-1] * 5
	
	state = 0
	
	runFile = open( runFileStr, 'r' )
	
	for line in runFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		#print( state, lineAr )
		if state == 0:
			if lineAr[0] == "Begin" and lineAr[1] == "converting":
				state = 1
			elif lineAr[0] == "Processed" and lineAr[1] == "reads:":
				outAr[0] += int( lineAr[2] )
			elif lineAr[0] == "Processed" and lineAr[1] == "bases:":
				outAr[1] += int(lineAr[2] )
			elif lineAr[0] == "Trimmed" and lineAr[1] == "reads:":
				outAr[2] += int( lineAr[2] )
				#print( lineAr[2] )
			elif lineAr[0] == "Quality-trimmed:":
				outAr[3] += int( lineAr[1] )
		elif state == 1 and lineAr[0] == "There":
			outAr[4] = int( lineAr[2] )
			state = 2
		elif state == 2 and lineAr[0] == "There":
			outAr[5] = int( lineAr[2] )
			outAr[6] = float( lineAr[5] )
			state = 3
		elif state == 3 and lineAr[2] == "p-value":
			outAr[7] = float( lineAr[8] )
			state = 4
		elif state == 4:
			#print( lineAr )
			outAr[8] = float(lineAr[4].replace("%",""))
			state = 5
			break
	runFile.close()
	return outAr

def readRunError( runErrorFileStr ):
	errorFile = open( runErrorFileStr, 'r' )
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
		elif line.startswith( '[mpileup]' ):
			break
	errorFile.close()
	return outAr

def readCorrentionOutput( corFileStr ):
	# lots-of cutoff and FDR
	outDict = {}
	nonConversion = -1
	
	state = 0
	
	corFile = open( corFileStr, 'r' )
	
	for line in corFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 5:
			continue
		if state == 0 and lineAr[1] == "non-conversion":
			nonConversion = float(lineAr[4].replace("%",""))
			state = 1
		elif state == 1 and lineAr[0] == "The" and lineAr[1] == 'closest':
			bases = lineAr[5]
			fdr = float( lineAr[-1] )
			outDict[bases] = fdr
		elif state == 1 and lineAr[0] == "Begin":
			break
	corFile.close()
	return nonConversion, outDict

def writeOutput( totalDict, cov, outFileStr ):
	# totalDict[n]=(run_array,correction_dict)
	outFile = open( outFileStr, 'w' )
	nucAr = ['CAA', 'CAC', 'CAT', 'CAN', 'CCA', 'CCC', 'CCT', 'CCN', 'CTA', 'CTC', 'CTT', 'CTN', 'CNA', 'CNC', 'CNT', 'CNN', 'CAG', 'CCG', 'CTG', 'CNG', 'CGA', 'CGC', 'CGG','CGT', 'CGN']
	headerStr = 'sample,processed reads,processed bases,adapter trimmed reads,quality trimmed bases,percent forward mapped,percent reverse mapped,uniquely mapped reads,non-clonal reads,percent reads remaining, genome pvalue,non-conversion percent'
	if cov != -1:
		headerStr += ',coverage'
	for i in nucAr:
		headerStr += ',FDR({:s})'.format(i)
	outFile.write( headerStr + '\n' )
	#print( totalDict['line1-G3-rep1'])
	#print( sorted(totalDict.keys()) )
	for sample in sorted(totalDict.keys()):
		tup = totalDict.get( sample )
		#try:
		ar = [ '{:d}'.format(x) for x in tup[0][0:4] ]
		ar += [ '{:.2f}'.format(x) for x in tup[0][4:6] ]
		ar += [ '{:d}'.format(x) for x in tup[0][6:8]]
		ar += [ '{:.5f}'.format(x) for x in tup[0][8:] ]
		if cov != -1:
			# mappedBases = uniquely mapped reads * floor( processBases / processReads )
			mappedBases = tup[0][6] * math.floor( float(tup[0][1]) / float(tup[0][0]) )
			ar += [ '{:.3f}'.format( mappedBases/cov) ]
		fdrAr = [ tup[1].get(x) for x in nucAr ]
		fdrAr = [('None' if x == None else '{:.5f}'.format(x)) for x in fdrAr ]
		
		outStr = sample+','+','.join(ar) + ',' + ','.join(fdrAr)+'\n'
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
	startInd = 0
	
	for i in range(min(3,len(argv)-2)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			fastaFileStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-u=' ):
			try:
				underScores = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of underscores must be an integer' )
				exit()
	# end for 
	pathRun = argv[startInd]
	pathCorrection = argv[startInd+1]
	
	processInputs( pathRun, pathCorrection, fastaFileStr, outID, underScores )

def printHelp():
	print ("Usage: python parse_pipeline_output.py [-f=fasta_file | fasta_index] [-o=out_id] [-u=underscore] <path_to_run_output> <path_to_correction_output>")
	print( 'Required: ')
	print( 'path_to_run_output\tpath to *.py.o* and *.py.e* files' )
	print( 'path_to_correction_output\tpath to *.sh.o* files' )
	print( 'Optional: ' )
	print( '-f=fasta\tfasta genome file or fasta index[default none]\n\t\tused to compute coverage [' )
	print( '-u=underscore_count\tunderscores is sample id name [default 0]' )
	print( '-o=out_id\toutput identifier' )
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
