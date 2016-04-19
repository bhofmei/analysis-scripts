import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 parse_rna_output_type2.py [-u=num_underscores_keep] [-o=out_prefix] <path_to_output>
# to be used without tophat alignment information
# this is when doing bowtie2 then cufflinks directly

UNDERSCORE=1

def processInputs( dirPath, outID, underScore ):
	totalDict = {}
	# get files
	outFiles = glob.glob( dirPath + '*.sh.e*' )
	
	print( 'Reading files for...' )
	for file in outFiles:
		nm = getSampleName( file, underScore )
		print( '-', nm )
		if totalDict.get(nm) == None:
			tempAr = readFile( file )
			totalDict[nm]=tempAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(nm, file) )
	
	if outID == None:
		outFileStr = 'rna_stats_out.csv'
	else:
		outFileStr = 'rna_{:s}_out.csv'.format( outID )
	writeOutput( outFileStr, totalDict )
	
			
def getSampleName( fileStr, underScore ):
	# file base name
	baseName = os.path.basename( fileStr )
	fInd = 0
	for i in range( underScore ):
		tInd = baseName[fInd+1:].find( '_' )
		if tInd != -1:
			fInd += tInd + 1
	if fInd ==0:
		tInd = baseName.find( '.' )
		return baseName[:tInd]
	return baseName[:fInd]

def readFile( fileStr ):
	
	# (0) input reads (1) surviving reads (2) dropped reads (3) aligned 0 times
	# (4) aligned 1 time (5) aligned >1 time (6) overall mapping rate 
	# (7) processed loci
	
	outAr = [-1] * 8
	state = 0
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
		elif state == 1 and lineAr[2] == "aligned":
			# aligned 0 times
			outAr[3] = int( lineAr[0] )
			state = 2
		elif state == 2 and lineAr[2] == "aligned":
			# aligned 1 time
			outAr[4] = int( lineAr[0] )
			state = 3
		elif state == 3 and lineAr[2] == "aligned":
			# aligned > 1 time
			outAr[5] = int( lineAr[0] )
			state = 4
		elif state == 4 and lineAr[2] == "alignment":
			outAr[6] = float( lineAr[0].replace("%","") )
			state = 5
		elif state == 5 and lineAr[1] == "Processed":
			outAr[7] = int( lineAr[2] )
			break
	# end for
	inFile.close()
	return outAr

def writeOutput( outFileStr, totalDict ):
	
	outFile = open( outFileStr, 'w' )
	headerStr = 'sample_name,input_reads,surviving_reads,dropped_reads,aligned_0_time,aligned_1_time,aligned_>1_times,overall_mapping_rate,processed_loci\n'
	outFile.write( headerStr )
	
	for sample in sorted( totalDict.keys() ):
		ar = totalDict[sample]
		outStr = sample
		for i in range(len(ar)):
			if i == 6:
				outStr += ',{:f}%'.format( ar[i] )
			else:
				outStr += ',{:d}'.format( ar[i] )
		outStr += '\n' 
		outFile.write( outStr )
	outFile.close()
		
def parseInputs( argv ):
	underScore = UNDERSCORE
	outID = None
	startInd = 0
	
	for i in range(min(2,len(argv)-1)):
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
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid argument'.format(argv[i]) )
			exit()
	# end for
	dirPath = argv[startInd]
	if dirPath.endswith( '/' ) == False:
		dirPath += '/'
	if os.path.isdir( dirPath ) == False:
		print( 'ERROR: {:s} is not a valid directory'.format( dirPath ) )
		exit()
	processInputs( dirPath, outID, underScore )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3.4 parse_rna_output_type2.py [-u=num_underscores_keep] [-o=out_id] <path_to_output>")
	else:
		parseInputs( sys.argv[1:] )
