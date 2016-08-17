import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python filter_pos_gbm_cds.p <pos_file> <gbm_file> <gff_file>

def processInputs(posFileStr, gbmFileStr, gffFileStr):
	info = '#from_script: filter_pos_gbm_cds.py; pos_file: {:s}; gbm_file: {:s}; gff_file: {:s}'.format( os.path.basename( posFileStr ), os.path.basename( gbmFileStr ), os.path.basename( gffFileStr ) )
	print( 'Position file:', os.path.basename( posFileStr ))
	print( 'gBM file:', os.path.basename( gbmFileStr ))
	print( 'GFF file:', os.path.basename( gffFileStr ) )
	base = os.path.basename( posFileStr )
	rInd = base.rfind( '.' )
	outFileStr = base[:rInd] + '_gbm' + base[rInd:]
	print( ' Reading gBM list' )
	gbmAr = readGBM( gbmFileStr )
	#print( gbmAr )
	print( ' Getting gBM CDS coordinates' )
	gbmDict = readGFF( gffFileStr, gbmAr )
	print( ' Filtering position file' )
	countB, countA = writeOutput( posFileStr, outFileStr, gbmDict, info )
	print( ' Original position count:', countB )
	print( ' gBM position count:', countA )
	print( ' Output written to', outFileStr )

def readGBM( gbmFileStr ):
	
	gbmFile = open( gbmFileStr, 'r' )
	gbmAr = []
	
	for line in gbmFile:
		if line.startswith( '#' ):
			continue
		line = line.rstrip()
		geneName = line
		gbmAr +=[ geneName ]
	gbmFile.close()
	
	return gbmAr

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

def readGFF( gffFileStr, gbmAr ):
	
	gffFile = open( gffFileStr, 'r' )
	gffDict = {}
	
	for line in gffFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		chrm = lineAr[0]
		strand = lineAr[6]
		if gffDict.get( chrm ) == None:
			gffDict[chrm] = []
		
		if lineAr[2] == "CDS":
			name = getGeneNameCDS( lineAr[8] )
			#print( '-', name )
			if name in gbmAr:
				for i in range(start, end+1 ):
					bisect.insort( gffDict[chrm], i )
		# end if CDS
	# end for line
	
	gffFile.close()
	return gffDict

def getGeneNameCDS (notesStr):
	search = "Parent="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]
	
def writeOutput( posFileStr, outFileStr, gbmDict, info ):
	countB = 0
	countA = 0
	inFile = open( posFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\nchrm\tpos\n' )
	for line in inFile:
		if line.startswith( '#' ) or line == '/n' or line.startswith( 'chrm' ):
			continue
		countB += 1
		lineAr = line.rstrip().split( '\t' )
		chrm = lineAr[0]
		pos = int( lineAr[1] )
		posAr = gbmDict.get(chrm)
		if posAr == None:
			continue
		isIn = bisectIndex( posAr, pos )
		if isIn != None:
			countA += 1
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()
	return countB, countA	

def parseInputs( argv ):
	posFileStr = argv[0]
	gbmFileStr = argv[1]
	gffFileStr = argv[2]
	processInputs(posFileStr, gbmFileStr, gffFileStr)


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print ("Usage: python filter_pos_gbm_cds.p <pos_file> <gbm_file> <gff_file>")
	else:
		parseInputs( sys.argv[1:] )
