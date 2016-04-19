import sys, glob

# Usage: python3.4 parse_rna_output.py <path_to_output> 

''' Assumes 1 file per sample where the first part of that file is the .sh.e output of the trimmomatic-tophat-cufflinks script and second part is align_summary.txt
'''

def readFile( inFileStr ):

	inFile = open( inFileStr, 'r' )
	# (0) input reads (1) surviving reads (2) dropped reads (3) processed loci
	# (4) aligned 1 time (5) aligned >1 time (6) aligned >20 times
	# (7) overall mapping rate
	outAr = [-1]*8
	
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
		elif state == 1 and lineAr[1] == "Processed":
			state = 2
		elif state ==1 and lineAr[0] == "Processed" and lineAr[2] == "loci.":
			outAr[3] = int( lineAr[1] )
			state = 3
		elif state == 2 and lineAr[1] == "Processed":
			outAr[3] = int( lineAr[2] )
			state = 3
		elif state == 3 and lineAr[0] == "Mapped":
			outAr[4] = int( lineAr[2] )
			state = 4
		elif state == 4:
			outAr[5] = int( lineAr[2] )
			#print( lineAr )
			try:
				outAr[6] =  int( lineAr[7].replace("(","") )
			except ValueError:
				outAr[6] =  int( lineAr[8].replace("(","") )
			state = 5
		elif state == 5:
			outAr[7] = float( lineAr[0].replace("%","") )
			break
	
	inFile.close()
	return outAr

def writeOutput( outFileStr, totalDict ):
	outFile = open( outFileStr, 'w' )
	headerStr = 'sample name,input reads,surviving reads,dropped reads,processed loci,aligned 1 time,aligned >1 time,aligned >20 times,overall mapping rate\n'
	outFile.write( headerStr )
	
	for sample in sorted(totalDict.keys()):
		ar = totalDict[sample]
		outStr = sample
		for i in range( len(ar) ):
			if i == 7:
				outStr += ',{:.1f}'.format( ar[i] )
			else:
				outStr += ',{:d}'.format( ar[i] )
		outStr += '\n'
		outFile.write( outStr )
	
	outFile.close()

def getSampleName( fileStr ):
	lind = fileStr.rindex( '/' )
	fileName = fileStr[lind+1:]
	return fileName.replace("run_info_","").replace(".txt","")

def mainFunction( pathOutput ):
	totalDict = {}
	outFiles = glob.glob(pathOutput+'*.txt')
	#print( outFiles )
	
	print( 'Reading files for...')
	for file in outFiles:
		n = getSampleName( file )
		print( '-',n)
		if totalDict.get(n) == None:
			tempAr = readFile( file )
			#print( tempAr )
			totalDict[n]=tempAr
		else:
			print( 'ERROR: duplicate run files for {:s}, skipping {:s}'.format(n, file) )
	print( 'Writing output...')
	writeOutput('rna_output.csv', totalDict )
	print( 'Done.' )

if __name__ == "__main__":
	if len(sys.argv) != 2 :
		print ("Usage: python3.4 parse_rna_output.py <path_to_output> ")
	else:
		pathOutput = sys.argv[1]
		mainFunction( pathOutput )
