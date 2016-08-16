import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3 parse_rna_output_type3.py [-u=num_underscores_keep] [-o=out_id] [-t=output_type] <path_to_output> 

# use for different rna-seq output types
# 1 - trimmomatic, bowtie2 SE (tbs)
# 2 - trimmomatic, tophat SE (tts)
## 3 - trimmomatic, bowtie2 PE (tbp)
## 4 - trimmomatic, tophat PE (ttp)
# 5 - trimmomatic, bowtie2, cufflinks SE (tbcs)
# 6 - trimmomatic, tophat, cufflinks SE [default] (ttcs)
## 7 - trimmomatic, bowtie2, cufflinks PE (tbcp)
# 8 - trimmomatic, tophat, cufflinks PE (tbcp)

UNDERSCORE=1
OUTTYPE=6

def processInputs( dirPath, outID, outType, underScore ):
	totalDict = {}
	outFiles = glob.glob( dirPath + '*.txt' )
	if len( outFiles ) == 0 :
		outFiles = glob.glob( dirPath + '*.sh.e*' )
		if len( outFiles ) == 0:
			print( 'ERROR: no files in directory or files incorrectly named' )
			exit()
	# loop through files
	print( 'Reading files' )
	for file in outFiles:
		nm = getSampleName( file, underScore )
		print( '-', nm )
		# read file based on type
		if totalDict.get(nm) == None:
			tmpAr = readFile( file, outType )
			totalDict[nm] = tmpAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(nm, file) )
	
	# output file name
	if outID == None:
		outFileStr = 'rna_stats_out.csv'
	else:
		outFileStr = 'rna_stats_{:s}.csv'.format( outID )
	print( 'Writing output to {:s}'.format( outFileStr ) )
	writeOutput( outFileStr, totalDict, outType )

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

def readFile( fileStr, outType ):
	if outType == 1: # TBS
		return readFile_tbp( fileStr )
	elif outType == 2: #TTS
		return readFile_tts( fileStr )
	elif outType == 3: # TBP - bad
		print( 'ERROR: output type 3 not currently supported' )
		exit()
	elif outType == 4: # TTP - bad
		print( 'ERROR: output type 4 not currently supported' )
		exit()
	elif outType == 5: # TBCS
		return readFile_tbcs( fileStr )
	elif outType == 6: #TTCS
		return readFile_ttcs( fileStr )
	elif outType == 7: # TBCP - bad
		print( 'ERROR: output type 7 not currently supported' )
		exit()
	elif outType == 8: # TTCP 
		return readFile_ttcp( fileStr )

def readFile_tbs( fileStr ):
	'''	TRIMMOMATIC, BOWTIE2 SE
	'''
	# (0) input reads (1) surviving reads (2) dropped reads 
	# (3) aligned 0 times (4) aligned 1 time (5) aligned >1 time 
	# (6) overall mapped (7) overall mapping rate (8) mulimapping rate
	
	outAr = [-1] * 9
	state = 0
	inFile = open( fileStr, 'r' )
	btInput = 0
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		
		if len(lineAr) < 3:
			continue
		if state == 0 and lineAr[0] == "Input":
			# Input Reads: 30868495 Surviving: 30835023 (99.89%) Dropped: 33472 (0.11%)
			outAr[0] = int( lineAr[2] )	# input reads
			outAr[1] = int( lineAr[4] )	# surviving reads
			outAr[2] = int( lineAr[7] )	# dropped reads
			state = 1
		if state == 1 and lineAr[1] == "reads;":
			# 60876989 reads; of these:
			btInput = int( lineAr[0] )
			state = 2
		elif state == 2 and lineAr[2] == "aligned":
			# aligned 0 times
			outAr[3] = int( lineAr[0] )
			state = 3
		elif state == 3 and lineAr[2] == "aligned":
			# aligned 1 time
			outAr[4] = int( lineAr[0] )
			state = 4
		elif state == 4 and lineAr[2] == "aligned":
			# aligned > 1 time
			outAr[5] = int( lineAr[0] )
			break	
	# end for line
	outAr[6] = outAr[4] + outAr[5] # overall mapped = aligned 1 + aligned > 1
	outAr[7] = float( outAr[6] ) / btInput * 100 # mapping rate = overall mapped / input
	outAr[8] = ( -1 if outAr[6] == 0 else (float( outAr[5] ) / outAr[6]) * 100) # multimapping = aligned > 1 / overall mapped
	inFile.close()
	return outAr
	
def readFile_tts( fileStr ):
	''' TRIMMOMATIC, TOPHAT SE
	'''
	# (0) input reads (1) surviving reads (2) dropped reads
	# (3) aligned 0 times (4) aligned 1 time (5) aligned >1 time 
	# (6) aligned >20 times (7) overall mapped 
	# (8) overall mapping rate (9) multimapping rate (10) overall remaining
	
	inFile = open( fileStr, 'r' )
	outAr = [-1]*11
	tpInput = 0
	
	state = 0
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		
		# skip lines shorter than 3 words
		if len(lineAr) < 3:
			continue
		if state == 0 and lineAr[0] == "Input":
			# Input Reads: 42200104 Surviving: 42125263 (99.82%) Dropped: 74841 (0.18%)
			outAr[0] = int( lineAr[2] )	# input reads
			outAr[1] = int( lineAr[4] )	# surviving reads
			outAr[2] = int( lineAr[7] )	# dropped reads
			state = 1
		elif state == 1 and lineAr[0] == "Input":
			#Input     :  42125263
			tpInput = int( lineAr[2] )
			state = 2
		elif state == 2:
			# Mapped   :  38836800 (92.2% of input)
			outAr[7] = int( lineAr[2] )
			state = 3
		elif state == 3:
			# of these:  33853373 (87.2%) have multiple alignments (3957 have >20)
			outAr[5] = int( lineAr[2] )
			try:
				outAr[6] = int( lineAr[7].replace('(','') )
			except ValueError:
				outAr[6] = int( lineAr[8].replace('(','') )
			break
	# end for line
	outAr[3] = 	tpInput - outAr[7] # aligned 0 = input - total mapped(7)
	outAr[4] = 	outAr[7] - outAr[5] # aligned 1 = total mapped(7) - aligned >1(5)
	outAr[8] = float(outAr[7]) / tpInput * 100 # overall mapping = total mapped (7) / input
	outAr[9] = float(outAr[5]) / outAr[7] * 100 # multimapping = >1 (5) / total mapped (7)
	outAr[10] = float(outAr[7]) / outAr[0] * 100 # remaining = total mapped (7) / input (0)
	inFile.close()
	return outAr

def readFile_tbcs( fileStr ):
	''' TRIMMOMATIC, BOWTIE2, CUFFLINKS SE
	'''
	# (0) input reads (1) surviving reads (2) dropped reads 
	# (3) aligned 0 times (4) aligned 1 time (5) aligned >1 time 
	# (6) overall mapped (7) overall mapping rate (8) mulimapping rate
	# (9) processed loci
	
	outAr = [-1] * 10
	state = 0
	btInput = 0
	inFile = open( fileStr, 'r' )
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		
		if len(lineAr) < 3:
			continue
		if state == 0 and lineAr[0] == "Input":
			# Input Reads: 30868495 Surviving: 30835023 (99.89%) Dropped: 33472 (0.11%)
			outAr[0] = int( lineAr[2] )	# input reads
			outAr[1] = int( lineAr[4] )	# surviving reads
			outAr[2] = int( lineAr[7] )	# dropped reads
			state = 1
		elif state == 1 and lineAr[1] == "reads;":
			btInput = int( lineAr[0] )
			state = 2
		elif state == 2 and lineAr[2] == "aligned":
			# aligned 0 times
			outAr[3] = int( lineAr[0] )
			state = 3
		elif state == 3 and lineAr[2] == "aligned":
			# aligned 1 time
			outAr[4] = int( lineAr[0] )
			state = 4
		elif state == 4 and lineAr[2] == "aligned":
			# aligned > 1 time
			outAr[5] = int( lineAr[0] )
			state = 4
		elif state == 4 and lineAr[1] == "Processed":
			outAr[9] = int( lineAr[2] )
			break
	# end for line
	outAr[6] = outAr[4] + outAr[5] # overall mapped = aligned 1 + aligned > 1
	outAr[7] = float( outAr[6] ) / btInput * 100 # mapping rate = overall mapped / input
	outAr[8] = float( outAr[5] ) / outAr[6] # multimapping = aligned > 1 / overall mapped
	
	inFile.close()
	return outAr

def readFile_ttcs( fileStr ):
	''' TRIMMOMATIC, TOPHAT, CUFFLINKS SE
	'''
	inFile = open( fileStr, 'r' )
	# (0) input reads (1) surviving reads (2) dropped reads
	# (3) aligned 0 times (4) aligned 1 time (5) aligned >1 time 
	# (6) aligned >20 times (7) overall mapped (8) mapping rate 
	# (9) multimapping rate (10) percent remaining (11) processed loci
	
	
	outAr = [-1]*12
	tpInput = 0
	state = 0
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		if state == 0 and lineAr[0] == "Input":
		# Input Reads: 42200104 Surviving: 42125263 (99.82%) Dropped: 74841 (0.18%)
			outAr[0] = int( lineAr[2] ) # input reads 
			outAr[1] = int( lineAr[4] ) # surviving reads
			outAr[2] = int( lineAr[7] ) # discarded reads
			state = 1
		elif state == 1 and lineAr[1] == "Processed":
			state = 2
		elif state ==1 and lineAr[0] == "Processed" and lineAr[2] == "loci.":
			outAr[11] = int( lineAr[1] )
			state = 3
		elif state == 2 and lineAr[1] == "Processed":
			outAr[11] = int( lineAr[2] )
			state = 3
		elif state == 3 and lineAr[0] == "Input":
			#Input     :  42125263
			tpInput = int( lineAr[2] )
			state = 4
		elif state == 4:
			# Mapped   :  38836800 (92.2% of input)
			outAr[7] = int( lineAr[2] )
			state = 5
		elif state == 5:
			# of these:  33853373 (87.2%) have multiple alignments (3957 have >20)
			outAr[5] = int( lineAr[2] )
			try:
				outAr[6] =  int( lineAr[7].replace("(","") )
			except ValueError:
				outAr[6] =  int( lineAr[8].replace("(","") )
			break
	# end for line
	outAr[3] = 	tpInput - outAr[7] # aligned 0 = input - total mapped(7)
	outAr[4] = 	outAr[7] - outAr[5] # aligned 1 = total mapped(7) - aligned >1(5)
	outAr[8] = float(outAr[7]) / tpInput * 100 # overall mapping = total mapped (7) / input
	outAr[9] = float(outAr[5]) / outAr[7] * 100 # multimapping = >1 (5) / total mapped (7)
	outAr[10] = float(outAr[7]) / outAr[0] * 100 # percent remaining = total mapped (7) / input (0)
	inFile.close()
	return outAr
	
def readFile_ttcp( fileStr ):
	''' TRIMMOMATIC, TOPHAT, CUFFLINKS PE
	'''
	inFile = open( fileStr, 'r' )
	# (0) input reads (1) surviving reads (2) dropped reads 
	# (3) pairs surviving (4) forward-only surviving (5) reverse-only surviving
	# (6) left overall mapped (7) left aligned 0 times 
	# (8) left aligned >20 times (9) right overall mapped
	# (10) right aligned 0 times (11) right aligned >20 times 
	# (12) overall mapped (13) overall mapping rate 
	# (14) aligned pairs (15) multiple alignment pairs (16) discordant pairs
	# (17) concordant pair rate (18) processed loci
	
	
	outAr = [-1]*19
	tpInput = [0,0]
	state = 0
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 2:
			continue
		if state == 0 and lineAr[0] == "Input":
		# Input Read Pairs: 40014844 Both Surviving: 39993127 (99.95%) Forward Only Surviving: 10936 (0.03%) Reverse Only Surviving: 1238 (0.00%) Dropped: 9543 (0.02%)
			outAr[0] = int( lineAr[3] ) # input reads 
			outAr[3] = int( lineAr[6] ) # pairs surviving reads
			outAr[4] = int( lineAr[11] ) # forward surviving reads
			outAr[5] = int( lineAr[16] ) # reverse surviving reads
			outAr[2] = int( lineAr[19] ) # discarded reads
			state = 1
		# error state for trimmomatic
		elif state == 0 and lineAr[0] == 'Exception':
			state = 1
		elif state == 1 and lineAr[1] == "Processed":
			state = 2
		elif state == 1 and lineAr[0] == "Processed" and lineAr[2] == "loci.":
			outAr[18] = int( lineAr[1] )
			state = 3
		elif state == 2 and lineAr[1] == "Processed":
			outAr[18] = int( lineAr[2] )
			state = 3
		elif state == 3 and lineAr[0] == 'Left':
			state = 4
		elif state == 3 and lineAr[0] == 'Right':
			state = 5
		elif (state == 4 or state == 5 ) and lineAr[0] == "Input":
			#Input     :  42125263
			q = state - 4
			tpInput[q] = int( lineAr[2] )
		elif (state == 4 or state == 5 ) and lineAr[0] == 'Mapped':
			# Mapped   :  38836800 (92.2% of input)
			# (6) left (9)
			q = (6 if state == 4 else 9 )
			outAr[q] = int( lineAr[2] )
		elif (state == 4 or state == 5 ):
			#  of these:  11019027 (28.9%) have multiple alignments (2571446 have >20)
			# (8) left (11) right
			q = (8 if state == 4 else 11 )
			try:
				outAr[q] =  int( lineAr[7].replace("(","") )
			except ValueError:
				outAr[q] =  int( lineAr[8].replace("(","") )
			state = 3
		elif state == 3 and lineAr[1] == 'overall':
			# 94.1% overall read mapping rate.
			outAr[13] = float( lineAr[0].replace('%','' ) )
			state = 6
		elif state == 6 and lineAr[0] == 'Aligned':
			# Aligned pairs:  36072458
			outAr[14] = int( lineAr[2] )
		elif state == 6 and lineAr[0] == 'of':
			# of these:  10583302 (29.3%) have multiple alignments
			outAr[15] = int( lineAr[2] )
			state = 7
		elif state == 7:
			# 2285255 ( 6.3%) are discordant alignments
			outAr[16] = int( lineAr[0] )
			state = 8
		elif state == 8 and lineAr[1] == 'concordant':
			# 84.5% concordant pair alignment rate.
			outAr[17] =  float( lineAr[0].replace('%','' ) )
			break
	# end for line
	outAr[1] = outAr[3] + outAr[4] + outAr[5] # surviving(1) = pairs(3) + forward($) + reverse(5) surviving
	outAr[7] = 	tpInput[0] - outAr[6] # left.aligned0 (7) = left.input - left.mapped(6)
	outAr[10] = tpInput[1] - outAr[9] # right.aligned0 (10) = right.input - right.mapped(9)
	outAr[12] = int(outAr[13] / 100 * ( tpInput[0] + tpInput[1] )) # overall mapped = mapping rate * (left.input + right.input)
	
	inFile.close()
	return outAr

def writeOutput( outFileStr, totalDict, outType ):
	btAr = ['sample_name', 'input_reads', 'surviving_reads', 'dropped_reads', 'aligned_0_times', 'aligned_1_time', 'aligned_>1_time', 'overall_mapped', 'overall_mapping_rate', 'multimapping_rate' ]
	tpAr = ['sample_name', 'input_reads', 'surviving_reads', 'dropped_reads', 'aligned_0_times', 'aligned_1_time', 'aligned_>1_time', 'aligned_>20_times', 'overall_mapped', 'overall_mapping_rate', 'multimapping_rate', 'percent_remaining' ]
	tpPeAr = ['sample_name','input_reads', 'surviving_reads', 'dropped_reads', 'pairs_surviving', 'forwand-only_surviving', 'reverse-only_surviving', 'left_overall_mapped', 'left_aligned_0_times', 'left_aligned_>20_times', 'right_overall_mapped', 'right_aligned_0_times', 'right_aligned_>20_times', 'overall_mapped', 'overall_mapping_rate', 'aligned_pairs', 'multiple_aligned_pairs', 'discordant_pairs', 'concordant_pair_rate' ]
	cfAr = [ 'processed_loci' ]
	# mod 2 means tophat, mod 4 == 0 or 1 means PE
	isBowtie = outType % 2
	isPE = outType % 4 # is PE if < 2
	if isBowtie:
		headerAr = btAr
	elif isPE < 2:	# ttp, ttcp
		headerAr = tpPeAr
	else: # tts, ttcs
		headerAr = tpAr
		
	# outtype > 5 means cufflinks
	if outType > 5:
		headerAr += cfAr
	
	outFile = open( outFileStr, 'w' )
	outFile.write( '{:s}\n'.format( ','.join( headerAr ) ) )
	
	ip = (8 if (outType % 2 == 0) else 7 )
	pefl = [ 13, 17 ]
	for sample in sorted( totalDict.keys() ):
		ar = totalDict[ sample ]
		#print( ar )
		outAr = [ sample ]
		if isPE > 1: # single end data
			outAr += [ '{:d}'.format( ar[i] ) for i in range( ip ) ]	# integers
			outAr += [ '{:.2f}'.format( ar[i] ) for i in range( ip, ip+2 ) ]
			if outType % 2 == 0:
				outAr += [ '{:.2f}'.format( ar[ip+2] ) ]
			if outType > 2:
				outAr += [ '{:d}'.format( ar[ip+3] ) ]
		else: # paired end
			outAr += [ ('{:.2f}'.format( ar[i] ) if i in pefl else '{:d}'.format( ar[i] ) ) for i in range(len(ar)) ]
		outFile.write( '{:s}\n'.format( ','.join( outAr ) ) )
	# end for sample	
	outFile.close()

def parseInputs( argv ):
	underScore = UNDERSCORE
	outID = None
	outType = OUTTYPE
	startInd = 0
	
	for i in range(min(3,len(argv))):
		if argv[i].startswith( '-u=' ):
			try:
				underScore = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of underscores to keep must be an integer' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			try:
				outType = int( argv[i][3:] )
				if outType < 1 or outType > 8:
					print( 'ERROR: not valid output type' )
					exit()
				startInd += 1
			except ValueError:
				print( 'ERROR: output type to keep must be an integer' )
				exit()
		elif argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid argument'.format(argv[i]) )
			exit()
	
	dirPath = argv[startInd]
	if dirPath.endswith( '/' ) == False:
		dirPath += '/'
	if os.path.isdir( dirPath ) == False:
		print( 'ERROR: {:s} is not a valid directory'.format( dirPath ) )
		exit()
	processInputs( dirPath, outID, outType, underScore )

def printHelp():
	print( 'Usage: python3 parse_rna_output.py [-u=num_underscores_keep] [-o=out_id] [-t=output_type] <path_to_output>' )
	print()
	print( "Required:\npath_to_output\tpath to directory with output files;\n\t\tall files must end with '.txt'\n\t\tor be of form '*.sh.e*'")
	print( "Optional:" )
	print( "-u=num_underscores\tnumber of underscores in sample name to keep;\n\t\t0 keeps all [default: 1]" )
	print( "-o=out_id\tidentifier used for output file [default: 'out']" )
	print( "-t=output_type\toutput type based on software steps used" )
	print( "\t\t1 - trimmomatic, bowtie2 SE" )
	print( "\t\t2 - trimmomatic, tophat SE" )
	print( "\t\t5 - trimmomatic, bowtie2, cufflinks SE" )
	print( "\t\t6 - trimmomatic, tophat, cufflinks SE [default]" )
	print( "\t\t7 - trimmomatic, tophat PE" )
	print( "\t\t8 - trimmomatic, tophat, cufflinks PE [default]" )
	print( "-h\tprint this help screen" )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
