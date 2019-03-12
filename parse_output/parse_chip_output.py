import sys, glob, os

# Usage: python3.4 parse_chip_output.py [-o=out_id] [-x=grep_pattern] [-u=underscore_keep] <path_to_output> 

UNDERSCORE=1

def mainFunction( pathOutput, outID, grepPattern, underScore ):
	totalDict = {}
	if pathOutput.endswith('/')==False:
		pathOutput += "/"
	outFiles = glob.glob(os.path.normpath(pathOutput+grepPattern))
	
	print( 'Reading files for...')
	for file in outFiles:
		n = getSampleName( file, underScore )
		print( '-',n)
		if totalDict.get(n) == None:
			tempAr = readFile( file )
			totalDict[n]=tempAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(n, file) )
	if outID == None:
		outFileStr = 'chip_output.csv'
	else:
		outFileStr = 'chip_output_' + outID + '.csv'
	print( 'Writing output to {:s}...'.format( outFileStr ))
	writeOutput( outFileStr, totalDict )
	print( 'Done.' )

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

def readFile( inFileStr ):

	inFile = open( inFileStr, 'r' )
	# (0) input reads (1) surviving reads (2) dropped reads (3) aligned 0 times 
	# (4) aligned 1 time (5) aligned >1 time (6) overall mapping rate
	# (7) mapped reads (8) clonal reads (9) percent duplication (10) remaining reads
	outAr = [-1]*11
	
	state = 0
	
	for line in inFile:
		line = line.rstrip().lstrip()
		lineAr = line.split()
		if len( lineAr ) < 3:
			continue
		if state == 0 and lineAr[0] == "Input":
			# input reads
			outAr[0] = int( lineAr[2] )
			# surviving reads
			outAr[1] = int( lineAr[4] )
			# discarded reads
			outAr[2] = int( lineAr[7] )
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
		elif state == 5 and lineAr[0] == "remove":
			state = 6
		### samtools remove
		elif state == 6 or (state == 5 and lineAr[0] == '[bam_rmdupse_core]'):
			try:
				outAr[7] = int( lineAr[3] )
				outAr[8] = int( lineAr[1] )
				outAr[9] = float( outAr[8] ) / float( outAr[7] ) * 100
				outAr[10] = outAr[7] - outAr[8]
			except ValueError:
				print( 'WARNING: no clonal duplicate removal information' )
			break
		### picard tools remove
		elif state == 5 and lineAr[0] == 'Unknown':
			outAr[7] = int( lineAr[2] )
			outAr[8] = int( lineAr[6] )
			outAr[9] = float( outAr[8] ) / float( outAr[7] ) * 100
			outAr[10] = outAr[7] - outAr[8]
			break
	
	inFile.close()
	return outAr

def writeOutput( outFileStr, totalDict ):
	outFile = open( outFileStr, 'w' )
	headerStr = 'sample name,input reads,surviving reads,dropped reads,aligned 0 times,aligned 1 time,aligned >1 time,overall mapping percent,mapped reads,clonal reads,percent duplication,remaining reads,percent remaining\n'
	outFile.write( headerStr )
	
	for sample in sorted(totalDict.keys()):
		ar = totalDict[sample]
		outStr = sample
		for i in range( len(ar) ):
			if i == 6 or i == 9:
				outStr += ',{:.2f}'.format( ar[i] )
			else:
				outStr += ',{:d}'.format( ar[i] )
		outStr += ',{:.2f}'.format( float(ar[10])/float(ar[0])*100 )
		outStr += '\n'
		outFile.write( outStr )
	
	outFile.close()

def parseInputs( argv ):
	startInd = 0
	outID = None
	grepPattern = '*.sh.*'
	underScore = UNDERSCORE
	
	for i in range(min(3,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-x=' ):
			grepPattern = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-u=' ):
			try:
				underScore = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of underscores to keep must be an integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	pathOutput = argv[startInd]
	mainFunction( pathOutput, outID, grepPattern, underScore )

def printHelp():
	print ("Usage: python3.4 parse_chip_output.py [-o=outID] [-x=grep_pattern] [-u=num_underscore_keep] <path_to_output> ")
	print( 'Required:' )
	print( 'path_to_output\tfolder with files for analysis' )
	print( 'Optional:' )
	print( '-o=out_id\tstring identifier for output file name [default none]' )
	print( '-x=grep_pattern\tpattern used to search for output files\n\t\t[default "*.sh.*]' )
	print( "-u=num_underscores\tnumber of underscores in sample name to keep;\n\t\t0 keeps all [default: 1]" )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
		
