import sys, math, glob, multiprocessing, subprocess, os

# Usage: python add_chrom_gene.py [-c=col_num_gene] <gff_file> <in_file1> [in_fileN]*

def processInputs( gffFileStr, inFileStrAr, colNum ):
	
	print( 'reading gff' )
	gffDict = readGFF( gffFileStr )
	
	for file in inFileStrAr:
		print( 'processing', file )
		processFile( file, colNum, gffDict )

def readGFF( gffFileStr ):
	gffDict = {}
	gffFile = open( gffFileStr, 'r' )
	
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) < 3:
			print( 'ERROR: line {:s}'.format( line ) )
			continue
		# (0) chrm (1) source (2) type (3) start (4) stop (5) ? (6) strand
		# (7) ? (8) notes
		if 'gene' in lineAr[2]:
			gNameStr = getGeneName( lineAr[8] )
			gNameAr = gNameStr.split( ',' )
			for g in gNameAr:
				gffDict[ g ] = lineAr[0]
	gffFile.close()
	return gffDict

def getGeneName( notesStr ):
	# if it has gene= use that as gene name otherwise use Name=
	search1 = "Alias="
	search2 = "Name="
	outName = ''
	index = notesStr.find(search1)
	adIndex = index + len(search1)
	if index != -1:
		endIndex = notesStr[adIndex:].find(';')
		if endIndex == -1:
			outName += notesStr[adIndex:] + ','
		else:
			outName += notesStr[adIndex:endIndex+adIndex] + ','
	
	index = notesStr.find(search2)
	adIndex = index + len(search2)
	if index != -1:
		endIndex = notesStr[adIndex:].find(';')
		if endIndex == -1:
			outName += notesStr[adIndex:] + ','
		else:
			outName += notesStr[adIndex:endIndex+adIndex] + ','
	return outName[:-1]
	
def processFile( fileStr, colNum, gffDict ):
	
	if fileStr.endswith( 'tsv' ):
		delim = '\t'
	elif fileStr.endswith( 'csv' ):
		delim = ','
	rInd = fileStr.rfind( '.' )
	outFileStr = fileStr[:rInd] + '_chrm' + fileStr[rInd:]
	Sherrie Hines
	readFile( fileStr, outFileStr, colNum, delim, gffDict )
	
def readFile( inFileStr, outFileStr, colNum, delim, gffDict ):
	
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in inFile:
		lineAr = line.rstrip().split( delim )
		# header type 1
		if line.startswith( '#' ):
			lineAr[0] = lineAr[0][1:]
			outFile.write( '#chrm{:s}{:s}\n'.format( delim, delim.join( lineAr ) ) )
		# header type 2
		elif line.startswith( 'gene' ):
			outFile.write( 'chrm{:s}{:s}\n'.format( delim, delim.join( lineAr ) ) )
		# gene
		else:
			gName = lineAr[ colNum ]
			chrm = gffDict.get( gName )
			outFile.write( '{:s}{:s}{:s}\n'.format( ( 'NA' if chrm == None else chrm ), delim, delim.join(lineAr) ) )
	# end for
	inFile.close()
	outFile.close()

def parseInputs( argv ):
	# Usage: python3.4 add_chrom_gene.py [-c=col_num_gene] <gff_file> <in_file1> [in_fileN]*
	
	colNum = 0
	startInd = 0
	
	for i in range(min(1, len(argv)-2)):
		if argv[i].startswith( '-c=' ):
			try:
				colNum = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: column number indicating gene must be integer' )
				exit()
	
	gffFileStr = argv[startInd]
	fileStrAr = []
	for j in range( startInd+1, len(argv) ):
		fileStrAr += [ argv[j] ]
	processInputs( gffFileStr, fileStrAr, colNum )
		
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3 add_chrom_gene.py [-c=col_num_gene] <gff_file> <in_file1> [in_fileN]*")
	else:
		parseInputs( sys.argv[1:] )
