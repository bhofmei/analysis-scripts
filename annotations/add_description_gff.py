import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 add_description_gff.py <description_file> <gff_file>

def processInputs( descFileStr, gffFileStr ):
	print( 'Reading description file...' )
	descDict = readDesc( descFileStr )
	gffName = os.path.basename( gffFileStr )
	rInd = gffName.rfind( '.' )
	gffBase = ( gffName if rInd == -1 else gffName[:rInd] )
	gffEnd = ( '.gff' if rInd == -1 else gffName[rInd:] )
	outFileStr = gffBase + "_descriptions" + gffEnd
	print( 'Reading GFF file.\nOutput written to {:s}'.format( outFileStr ) )
	readGFF( gffFileStr, descDict, outFileStr )

def readDesc( descFileStr ):
	outDict = {}
	descFile = open( descFileStr, 'r' )
	isHeader = True
	for line in descFile:
		if isHeader:
			isHeader = False
			continue
		lineAr = line.rstrip().split('\t')
		# (0) Model_name (1) Type (2) Short_description (3) Curator_summary	
		# (4) Computational_description
		if len(lineAr) < 3:
			# (0) model name (1) description
			geneName = lineAr[0]
			rInd = geneName.rfind( '.' )
			if rInd != -1:
				geneName = geneName[:rInd]
			desc = lineAr[1]
			if desc != "":
				outDict[ geneName ] = desc
		else:
			geneName = lineAr[0]
			rInd = geneName.rfind( '.' )
			if rInd != -1:
				geneName = geneName[:rInd]
			desc = lineAr[2]
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
	descFileStr = argv[0]
	gffFileStr = argv[1]
	processInputs( descFileStr, gffFileStr )

if __name__ == "__main__":
	if len(sys.argv) != 3 :
		print ("python3.4 add_description_gff.py <description_file> <gff_file>")
	else:
		parseInputs( sys.argv[1:] )
