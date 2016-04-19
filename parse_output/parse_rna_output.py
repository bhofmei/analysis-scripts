import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3 parse_rna_output_type3.py [-u=num_underscores_keep] [-o=out_id] [-t=output_type] <path_to_output> 

# use for different rna-seq output types
# 1 - trimmomatic, bowtie2
# 2 - trimmomatic, tophat
# 3 - trimmomatic, bowtie2, cufflinks
# 4 - trimmomatic, tophat, cufflinks [default]

UNDERSCORE=1
OUTTYPE=4

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
	if outType == 1:
		return readFile1( fileStr )
	elif outType == 2:
		return readFile2( fileStr )
	elif outType == 3:
		return readFile3( fileStr )
	elif outType == 4:
		return readFile4( fileStr )

def readFile1( fileStr ):
	'''	TRIMMOMATIC, BOWTIE2
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
	
def readFile2( fileStr ):
	''' TRIMMOMATIC, TOPHAT
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

def readFile3( fileStr ):
	''' TRIMMOMATIC, BOWTIE2, CUFFLINKS
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

def readFile4( fileStr ):
	''' TRIMMOMATIC, TOPHAT, CUFFLINKS
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

def writeOutput( outFileStr, totalDict, outType ):
	btAr = ['sample_name', 'input_reads', 'surviving_reads', 'dropped_reads', 'aligned_0_times', 'aligned_1_time', 'aligned_>1_time', 'overall_mapped', 'overall_mapping_rate', 'multimapping_rate' ]
	tpAr = ['sample_name', 'input_reads', 'surviving_reads', 'dropped_reads', 'aligned_0_times', 'aligned_1_time', 'aligned_>1_time', 'aligned_>20_times', 'overall_mapped', 'overall_mapping_rate', 'multimapping_rate', 'percent_remaining' ]
	cfAr = [ 'processed_loci' ]
	if outType % 2 == 0:
		headerAr = tpAr
	else:
		headerAr = btAr
	if outType > 2:
		headerAr += cfAr
	
	outFile = open( outFileStr, 'w' )
	outFile.write( '{:s}\n'.format( ','.join( headerAr ) ) )
	# loop through samples
	ip = (8 if (outType % 2 == 0) else 7 )
	for sample in sorted( totalDict.keys() ):
		ar = totalDict[ sample ]
		outAr = [ sample ]
		
		outAr += [ '{:d}'.format( ar[i] ) for i in range( ip ) ]	# integers
		outAr += [ '{:.2f}'.format( ar[i] ) for i in range( ip, ip+2 ) ]
		if outType % 2 == 0:
			outAr += [ '{:.2f}'.format( ar[ip+2] ) ]
		if outType > 2:
			outAr += [ '{:d}'.format( ar[ip+3] ) ]
		outFile.write( '{:s}\n'.format( ','.join( outAr ) ) )
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
				if outType < 1 or outType > 4:
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
	print( "\t\t1 - trimmomatic, bowtie2" )
	print( "\t\t2 - trimmomatic, tophat" )
	print( "\t\t3 - trimmomatic, bowtie2, cufflinks" )
	print( "\t\t4 - trimmomatic, tophat, cufflinks [default]" )
	print( "-h\tprint this help screen" )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
