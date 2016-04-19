import sys, os, glob, math

# Usage: python3.4 parse_pipeline_output.py <path_to_run_output> <path_to_correction_output>

UNDERSCORE=3

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
		elif state == 1 and lineAr[0] == "The":
			bases = lineAr[5]
			fdr = float( lineAr[-1] )
			outDict[bases] = fdr
		elif state == 1 and lineAr[0] == "Begin":
			break
	corFile.close()
	return nonConversion, outDict

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

def readIndex( fastaIndexStr ):
	baseCount = 0
	indexFile = open( fastaIndexStr, 'r' )
	for line in indexFile:
		lineAr = line.rstrip().split('\t')
		# (0) name (1) length (2) byte position (3) line length w/o \n
		# (4) line length w/ \n
		if lineAr[0] not in ['mitochondria','chloroplast']:
			baseCount += int( lineAr[1] )
	indexFile.close()
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
		try:
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
		except TypeError:
			print( sample )
			print( tup )
			print( ar )
			print( fdrAr )
			exit()
		outFile.write( outStr )
	outFile.close()

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
			

def mainFunction( pathRun, pathCorrection, fastaFileStr, fastaIndexStr ):
	underscore = UNDERSCORE
	totalDict = {}
	runFiles = glob.glob(pathRun+'*.py.o*')
	errorFiles = glob.glob(pathRun+'*.py.e*')
	corFiles = glob.glob(pathCorrection+'*.sh.*')
	
	cov = -1
	if fastaFileStr != None:
		cov = readFasta( fastaFileStr )
		print( 'Read', fastaFileStr )
	elif fastaIndexStr != None:
		cov = readIndex( fastaIndexStr )
		print( 'Read', fastaIndexStr )
	print( 'Reading run files for...')
	for run in runFiles:
		n = getSampleName( run, underscore )
		print( '-',n)
		if totalDict.get(n) == None:
			tempAr = readRunOutput( run )
			totalDict[n]=tempAr
			#print( tempAr )
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(n, run) )
	print( 'Reading error files for...' )
	for err in errorFiles:
		n = getSampleName( err, underscore )
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
		n = getSampleName( cor, underscore )
		print('-',n )
		if totalDict.get(n) == None:
			print( 'ERROR: {:s} has correction file but not run file'.format(n ))
		else:
			tAr = totalDict.get(n)
			nonConversion, tempDict = readCorrentionOutput( cor )
			#print( tempDict )
			#if tAr[-1] == -1:
				#tAr[-1] = nonConversion
			totalDict[n] = ( tAr, tempDict )
	#print( totalDict )
	# write output
	print( 'Writing output...')
	writeOutput( totalDict, cov, 'run_correction_output.csv' )
	print( 'Done.' )
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 parse_pipeline_output.py [fasta_file | fasta_index] <path_to_run_output> <path_to_correction_output>")
	else:
		if sys.argv[1].endswith('.fa') or sys.argv[1].endswith('.fasta'):
			fastaFileStr = sys.argv[1]
			fastaIndexStr = None
			pathRun = sys.argv[2]
			pathCorrection = sys.argv[3]
		elif sys.argv[1].endswith('.fai'):
			fastaFileStr = None
			fastaIndexStr = sys.argv[1]
			pathRun = sys.argv[2]
			pathCorrection = sys.argv[3]
		else:
			fastaFileStr = None
			fastaIndexStr = None
			pathRun = sys.argv[1]
			pathCorrection = sys.argv[2]
		mainFunction( pathRun, pathCorrection, fastaFileStr, fastaIndexStr )
