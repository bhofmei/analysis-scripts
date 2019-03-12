import sys, os, bisect

# Usage: python filter_pos_gbm.py [-h] [-q] [-cds] [-v] <pos_file> <gbm_file> <gff_file>

def processInputs( posFileStr, subsetFileStr, gffFileStr, isCDS, isOppo, isPrint ):

	info = '#from_script: filter_pos_gene_gbm.py; pos_file: {:s}; subset_file: {:s}; gff_file: {:s}; use_cds: {:s}; is_opposite: {:s}'.format( os.path.basename( posFileStr ), ('None' if subsetFileStr == None else os.path.basename( subsetFileStr ) ), os.path.basename( gffFileStr ), str(isCDS), str(isOppo) )
	if isPrint:
		print( 'Position file:', os.path.basename( posFileStr ))
		print( 'Gene subset file:', ('None' if subsetFileStr == None else os.path.basename( subsetFileStr )))
		print( 'GFF file:', os.path.basename( gffFileStr ) )
		print( 'Using features: {:s}'.format(( 'CDS' if isCDS else 'gene' )) )
		print( 'Getting opposite coordinates:', isOppo )
	base = os.path.basename( posFileStr )
	rInd = base.rfind( '.' )

	genetype, cdstype = getOutputType( subsetFileStr, isCDS, isOppo )

	outFileStr = base[:rInd] + '_' + genetype + '_' + cdstype + base[rInd:]

	if subsetFileStr != None:
		if isPrint:
			print( ' Reading gene subset list' )
		geneAr = readSubset( subsetFileStr )
	else:
		geneAr = []
	if isPrint:
		print( ' Getting {:s} {:s} coordinates'.format(genetype, cdstype) )
	geneDict = readGFF( gffFileStr, geneAr, isCDS )
	if isPrint:
		print( ' Filtering position file' )
	countB, countA = writeOutput( posFileStr, outFileStr, geneDict, info, isOppo )
	if isPrint:
		print( ' Original position count:', countB )
		print( ' output position count:', countA )
		print( ' Output written to', outFileStr )

def getOutputType( subsetFileStr, isCDS, isOppo ):
	genetype = ''
	cdstype = 'CDS' if isCDS else 'gene'
	if subsetFileStr == None:
		genetype = 'nongenes' if isOppo else 'genes'
	else:
		genetype = 'oppo-subset' if isOppo else 'subset'

	return genetype, cdstype

def readSubset( inFileStr ):

	inFile = open( inFileStr, 'r' )
	outAr = []

	for line in inFile:
		if line.startswith( '#' ):
			continue
		line = line.rstrip()
		geneName = line
		outAr +=[ geneName ]
	inFile.close()

	return outAr

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

def readGFF( gffFileStr, geneAr, isCDS ):

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
			if len(geneAr) == 0 or name in geneAr:
				for i in range(start, end+1 ):
					bisect.insort( gffDict[chrm], i )
		# end if CDS
		elif lineAr[2] == 'mRNA':
			name = getGeneName( lineAr[8] )
			#print( '-', name )
			if len(geneAr) == 0 or name in geneAr:
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

def writeOutput( posFileStr, outFileStr, geneDict, info, isOppo ):
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
		posAr = geneDict.get(chrm)
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
	isPrint = True
	startInd = 0

	for i in range(min(3,len(argv)-3)):
		if argv[i].lower() == '-cds':
			isCDS = True
			startInd += 1
		elif argv[i] == '-v':
			isOppo = True
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for i

	posFileStr = argv[startInd]
	subsetFileStr = argv[startInd+1]
	# check for none/na condition
	if subsetFileStr.lower() in ['none', 'na']:
		subsetFileStr = None
	gffFileStr = argv[startInd+2]
	processInputs(posFileStr, subsetFileStr, gffFileStr, isCDS, isOppo, isPrint)

def printHelp():
	print ("Usage:\tpython filter_pos_gene_subset.py [-h] [-q] [-cds] [-v] <pos_file>\n\t<gbm_file> <gff_file>")
	print()
	print('Required:')
	print('pos_file\t1-based position file, tab-delimited columns: chrm, start, end')
	print('subset_file\tfile with list of genes to subset, one gene per line\n\t\tuse "none" or "na" to use all genes')
	print('gff_file\tGFF formatted file with genes')
	print()
	print('Optional:')
	print( '-h\t\tprint help and exit' )
	print( '-q\t\tquiet; do not print progress' )
	print('-cds\t\tuse CDS annotation not gene')
	print('-v\t\tinclude coordinates opposite of what is specified')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
