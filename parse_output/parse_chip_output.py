import sys, glob, os

# Usage: python3.4 parse_chip_output.py <path_to_output> 

def readFile( inFileStr ):

	inFile = open( inFileStr, 'r' )
	# (0) input reads (1) surviving reads (2) dropped reads (3) aligned 0 times 
	# (4) aligned 1 time (5) aligned >1 time (6) overall mapping rate
	# (7) mapped reads (8) clonal reads (9) remaining reads
	outAr = [-1]*10
	
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
		elif state == 6:
			try:
				outAr[7] = int( lineAr[3] )
				outAr[8] = int( lineAr[1] )
				outAr[9] = outAr[7] - outAr[8]
			except ValueError:
				print( 'WARNING: no clonal duplicate removal information' )
			break
	
	inFile.close()
	return outAr

def writeOutput( outFileStr, totalDict ):
	outFile = open( outFileStr, 'w' )
	headerStr = 'sample name,input reads,surviving reads,dropped reads,aligned 0 times,aligned 1 time,aligned >1 time,overall mapping percent,mapped reads,clonal reads,remaining reads,percent remaining\n'
	outFile.write( headerStr )
	
	for sample in sorted(totalDict.keys()):
		ar = totalDict[sample]
		outStr = sample
		for i in range( len(ar) ):
			if i == 6:
				outStr += ',{:.2f}'.format( ar[i] )
			else:
				outStr += ',{:d}'.format( ar[i] )
		outStr += ',{:.2f}'.format( float(ar[9])/float(ar[0])*100 )
		outStr += '\n'
		outFile.write( outStr )
	
	outFile.close()

def getSampleNameHelper( fileStr ):
	lind = fileStr.rindex( '/' )
	fileName = fileStr[lind+1:]
	rind = fileName.index('_')
	rind2 = fileName[rind+1:].index('_')
	rind3 = fileName[rind+rind2+2:].index('_')
	# return fileName[:rind2+rind+1]
	return fileName[:rind3+rind2+rind+2]

def getSampleName( fileStr ):
	fileName = os.path.basename( fileStr )
	# search for run
	runInd = fileName.find( 'run' )
	if runInd == -1:
		return getSampleNameHelper( fileStr )
	uInd = fileName.find( '_', runInd )
	return fileName[:uInd]
	

def mainFunction( pathOutput, outID ):
	totalDict = {}
	if pathOutput.endswith('/')==False:
		pathOutput += "/"
	outFiles = glob.glob(pathOutput+'*.sh.*')
	
	print( 'Reading files for...')
	for file in outFiles:
		n = getSampleName( file )
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

def parseInputs( argv ):
	startInd = 0
	outID = None
	for i in range(min(1,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
	pathOutput = argv[startInd]
	mainFunction( pathOutput, outID )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3.4 parse_chip_output.py [-o=outID] <path_to_output> ")
	else:
		parseInputs( sys.argv[1:] )
		
