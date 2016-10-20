import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python filter_pos_gbm.p [-cds] [-v] <pos_file> <gbm_file> <gff_file>

def processInputs(posFileStr, gbmFileStr, gffFileStr, isCDS, isOppo):
	info = '#from_script: filter_pos_gbm.py; pos_file: {:s}; gbm_file: {:s}; gff_file: {:s}; use_cds: {:s}; is_opposite: {:s}'.format( os.path.basename( posFileStr ), os.path.basename( gbmFileStr ), os.path.basename( gffFileStr ), str(isCDS), str(isOppo) )
	print( 'Position file:', os.path.basename( posFileStr ))
	print( 'gBM file:', os.path.basename( gbmFileStr ))
	print( 'GFF file:', os.path.basename( gffFileStr ) )
	print( 'Using gBM features: {:s}'.format(( 'CDS' if isCDS else 'gene' )) )
	print( 'Getting opposite coordinates:', isOppo )
	base = os.path.basename( posFileStr )
	rInd = base.rfind( '.' )
	outFileStr = base[:rInd] + '_' + ('ngbm' if isOppo else 'gbm') + '_' + ('cds' if isCDS else 'gene') + base[rInd:]
	print( ' Reading gBM list' )
	gbmAr = readGBM( gbmFileStr )
	#print( gbmAr )
	print( ' Getting {:s} {:s} coordinates'.format(('ngBM' if isOppo else 'gBM'), ('CDS' if isCDS else 'gene')  ) )
	gbmDict = readGFF( gffFileStr, gbmAr, isCDS )
	print( ' Filtering position file' )
	countB, countA = writeOutput( posFileStr, outFileStr, gbmDict, info, isOppo )
	print( ' Original position count:', countB )
	print( ' output position count:', countA )
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

def readGFF( gffFileStr, gbmAr, isCDS ):
	
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
			
		if lineAr[2] == "CDS" and isCDS:
			name = getGeneNameCDS( lineAr[8] )
			#print( '-', name )
			if name in gbmAr:
				for i in range(start, end+1 ):
					bisect.insort( gffDict[chrm], i )
		# end if CDS
		elif lineAr[2] == 'mRNA':
			name = getGeneName( lineAr[8] )
			#print( '-', name )
			if name in gbmAr:
				for i in range(start, end+1 ):
					bisect.insort( gffDict[chrm], i )
		# end elif mRNA
	# end for line
	
	gffFile.close()
	return gffDict

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def getGeneNameCDS (notesStr):
	search = "Parent="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]
	
def writeOutput( posFileStr, outFileStr, gbmDict, info, isOppo ):
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
		if (isIn != None and not isOppo) or (isIn == None and isOppo):
			countA += 1
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()
	return countB, countA	

def parseInputs( argv ):

	isCDS = False
	isOppo = False
	startInd = 0
	
	for i in range(min(2,len(argv)-3)):
		if argv[i].lower() == '-cds':
			isCDS = True
			startInd += 1
		elif argv[i] == '-v':
			isOppo = True
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for i
	
	posFileStr = argv[startInd]
	gbmFileStr = argv[startInd+1]
	gffFileStr = argv[startInd+2]
	processInputs(posFileStr, gbmFileStr, gffFileStr, isCDS, isOppo)


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print ("Usage: python filter_pos_gbm.py [-cds] [-v] <pos_file> <gbm_file> <gff_file>")
	else:
		parseInputs( sys.argv[1:] )
