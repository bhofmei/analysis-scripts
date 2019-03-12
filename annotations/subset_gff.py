import sys, os, bisect

# Usage: python subset_gff.py [-v] [-c=col_index] [-o=out_id] <gene_list_file> <gff_file>

def processInputs(gffFileStr, geneFileStr, outID, colIndex, isOppo):
	
	# read gene file
	print( 'Reading gene list' )
	geneList = readGeneFile( geneFileStr, colIndex )
	print( 'Found', len(geneList), 'genes' )
	
	if outID == None:
		i = os.path.basename(gffFileStr).rfind('.')
		if i != -1:
			outID = os.path.basename(gffFileStr)[:i]
			outID += '_subset'
	outFileStr = outID + '.gff'
	
	print( 'Writing output to', outFileStr )
	readGFF( gffFileStr, outFileStr, geneList, isOppo )
	print( 'Done' )
	
def readGeneFile( geneFileStr, colIndex ):
	# this file should be
	outAr = []
	inFile = open( geneFileStr, 'r' )
	
	for line in inFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith( '#' ):
			continue
		elif len(lineAr) < (colIndex+1):
			continue
		# otherwise its a good gene. remove .# if necessary
		geneName = lineAr[colIndex]
		i = geneName.rfind('.')
		if i != -1:
			geneName = geneName[:i]
		#outAr += [geneName]
		bisect.insort_left( outAr, geneName )
	# end for line
	inFile.close()
	return outAr

def readGFF( gffFileStr, outFileStr, geneList, isOppo ):
	inFile = open( gffFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	isWrite = False
	
	for line in inFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) < 9:
			continue
			
		fType = lineAr[2]
		# come across a gene
		if fType == 'gene':
			# get gene name
			gName = getGeneName( lineAr[8], 'Name=' )
			# check if in geneList
			j = bisectIndex( geneList, gName )
			if isOppo:
				isWrite = (j == None)
			else:
				isWrite = (j != None)
		
		if isWrite:
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()
			
def getGeneName( notesStr, search ):

	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	n = notesStr[adIndex:]
	if endIndex != -1:
		n = notesStr[adIndex:endIndex+adIndex]
		
	pInd = n.rfind( '.' )
	if pInd != 1:
		return n[:pInd]
	else:
		return n

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None
			
def parseInputs( argv ):
	outID = None
	isOppo = False
	colIndex = 0
	startInd = 0
	
	for i in range(min(3, len(argv)-2)):
		if argv[i] == '-v':
			isOppo = True
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			try:
				colIndex = int(argv[i][3:])
			except ValueError:
				print( 'ERROR: column index is not an integer...using default 0' )
			startInd += 1
	# end for
	geneFileStr = argv[startInd]
	gffFileStr = argv[startInd+1]
	processInputs( gffFileStr, geneFileStr, outID, colIndex, isOppo )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python subset_gff.py [-v] [-c=col_index] [-o=out_id] <gene_list_file> <gff_file>\n-v\t\texclude genes in file\n-c=col_index\t0-based index for column which has gene name [default 0]\n-o=out_id\talternative output file name [default 'input'_subset.gff]")
	else:
		parseInputs( sys.argv[1:] )
