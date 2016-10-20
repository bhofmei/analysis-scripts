import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python add_description_gff.py [-c=col_nums] <description_file> <gff_file>

COL=[0,1]

def processInputs( descFileStr, gffFileStr, colNums ):
	print( 'Reading description file...' )
	descDict = readDesc( descFileStr, colNums )
	gffName = os.path.basename( gffFileStr )
	rInd = gffName.rfind( '.' )
	gffBase = ( gffName if rInd == -1 else gffName[:rInd] )
	gffEnd = ( '.gff' if rInd == -1 else gffName[rInd:] )
	outFileStr = gffBase + "_descriptions" + gffEnd
	print( 'Reading GFF file.\nOutput written to {:s}'.format( outFileStr ) )
	readGFF( gffFileStr, descDict, outFileStr )

def readDesc( descFileStr, colNums ):
	outDict = {}
	descFile = open( descFileStr, 'r' )
	isHeader = True
	colNumsMax = max(colNums)
	for line in descFile:
		if isHeader:
			isHeader = False
			continue
		lineAr = line.rstrip().split('\t')
		if colNumsMax > len(lineAr):
			continue
		geneName = lineAr[colNums[0]]
		rInd = geneName.rfind( '.' )
		if rInd != -1:
			geneName = geneName[:rInd]
		desc = lineAr[colNums[1]]
		if desc != "":
			outDict[ geneName ] = desc
	# end for line
	descFile.close()
	return outDict

def readGFF( gffFileStr, descDict, outFileStr ):
	
	outFile = open( outFileStr, 'w' )
	gffFile = open( gffFileStr, 'r' )
	
	for line in gffFile:
		# headers
		if line.startswith( '#' ):
			outFile.write( line )
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		if lineAr[2] == "gene":
			geneName = getGeneName( lineAr[8] )
			desc = descDict.get( geneName )
			if desc != None:
				fDesc = desc.replace( ' ', '_' )
				lineAr[8] += ";Description=" + fDesc
		# write to file
		outFile.write( '\t'.join( lineAr ) + '\n' )
	# end for line
	outFile.close()
	gffFile.close()

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]
	

def parseInputs( argv ):
	colNums = COL
	startInd = 0
	for i in range(min(1,len(argv)-2)):
		if argv[i].startswith( '-c=' ):
			tmp = argv[i][3:].split(',')
			if len(tmp) != 2:
				print( 'ERROR: must specify both gene and description columns' )
				exit()
			try:
				colNums = [ int(x) for x in tmp ]
			except ValueError:
				print( 'ERROR: column indexes must be integers' )
				exit()
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for i
	
	descFileStr = argv[startInd]
	gffFileStr = argv[startInd+1]
	processInputs( descFileStr, gffFileStr, colNums )

def printHelp():
	print ("python3 add_description_gff.py [-c=col_nums] <description_file> <gff_file>")
	print( 'Adds functional descriptions to gene features of GFF' )
	print( 'Required:' )
	print( 'description_file\ttab-delimited file with, at least, gene name\n\t\tand functional description' )
	print( 'gff_file\tGFF file to add descriptions to' )
	print( 'Optional:' )
	print( '-c=col_nums\t0-index column numbers of gene name and description, respectively, comma-separated [default 0,1]' )
	
if __name__ == "__main__":
	if len(sys.argv) != 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
