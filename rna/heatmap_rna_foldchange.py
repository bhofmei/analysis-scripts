import sys, math, glob, multiprocessing, subprocess, os
ADD_VAL=0.5
HASHVAL=17
PERCENTILE=.90

# Usage: python3.4 heatmap_rna_foldchange.py [-h=h3vo_gene_file] [-j=base_j_gene_file] [-t=threshold] <gene_file> <combined_fpkm_file> <base_name> <compare_name> [compare_name]*

def readGeneFile( geneFileStr ):
	geneFile = open( geneFileStr, 'r' )
	geneAr = []
	isHeader = True
	
	for line in geneFile:
		if isHeader:
			isHeader = False
			continue
		line = line.rstrip().lstrip()
		if line != '':
			geneAr += [ line ]
	geneFile.close()
	return geneAr

def readFPKMFile( fpkmFileStr, sampleNamesAr, geneAr ):
	
	fpkmFile = open( fpkmFileStr, 'r' )
	geneIndex = None
	getColumns = None
	# fpkmMatrix arraged [gene][sample] - gene is rows, samples in columns
	fpkmMatrix = [ [-1]*len(sampleNamesAr) for i in range(len( geneAr )) ]
	
	for line in fpkmFile:
		if geneIndex == None:
			getColumns, geneIndex = parseHeader( line.rstrip(), sampleNamesAr )
			#print( getColumns)
			continue
		lineAr = line.rstrip().split( '\t' )
		geneName = lineAr[ geneIndex ]
		if geneName in geneAr:
			gInd = geneAr.index( geneName )
			for j in range( len( getColumns ) ):
				fpkmMatrix[gInd][j] = float( lineAr[getColumns[j]] )
	# finished parsing
	fpkmFile.close()
	return fpkmMatrix

def computeFoldChange( fpkmMatrix):
	
	outMatrix = [ [-1] * ( len(fpkmMatrix[0])-1 ) for x in range( len(fpkmMatrix))]
	numList = []
	
	# loop through genes
	for i in range( len(fpkmMatrix) ):
		# loop through comparison samples
		for j in range( 1, len(fpkmMatrix[i] ) ):
			fc = ( fpkmMatrix[i][j] + ADD_VAL ) / ( fpkmMatrix[i][0] + ADD_VAL )
			outMatrix[i][j-1] = fc
			numList += [fc]
	return outMatrix, numList

def calculatePercentile( percentile, numList ):
	
	if percentile == 1:
		return max( numList )
	numList.sort()
	ind = math.ceil( percentile * len( numList ) - 1 )
	p = numList[ind]
	if p == 0:
		print( '***** not percentile corrected *****' )
		return numList[-1]
	return p

def correctFoldChange( fpkmMatrix, threshold ):
	
	for i in range(len(fpkmMatrix)):
		for j in range(len(fpkmMatrix[i])):
			if fpkmMatrix[i][j] > threshold:
				fpkmMatrix[i][j] = threshold
	return fpkmMatrix

def orderGenes( fcMatrix, geneAr ):
	
	orderAr = []
	#print( fcMatrix)
	rankMatrix = rankTransform( fcMatrix )
	# loop through genes
	for i in range(len(fcMatrix)):
		hashVal = 0
		# loop through columns
		for j in range(len(fcMatrix[i])):
			# use ranked value
			#hashVal += j * rankMatrix[i][j] * len(fcMatrix[i])
			hashVal += j * rankMatrix[i][j] * HASHVAL
		orderAr += [(hashVal,i)]
	#print( 'orderAr\n',orderAr )
	ordered = sorted( orderAr )
	#print( ordered )
	outMatrix = []
	outGenes = []
	for tup in ordered:
		outMatrix.append( fcMatrix[ tup[1] ] )
		outGenes.append( geneAr[tup[1]] )
	#print( outMatrix, outGenes)
	return outMatrix, outGenes

def rankTransform( inMatrix ):
	# rank by column
	tmpMatrix = [ [-1] * len(inMatrix) for x in range( len(inMatrix[0] ) ) ]
	outMatrix = [ [-1]*len(inMatrix[0]) for x in range(len(inMatrix)) ]
	# loop through rows
	for i in range(len(inMatrix) ):
		# loop through columns
		for j in range( len(inMatrix[i]) ):
			tmpMatrix[j][i] = (inMatrix[i][j],i)
	# tmpMatrix = [(change,index)]
	# sort rows of tmpMatrix
	for j in range(len(tmpMatrix)):
		row = tmpMatrix[j]
		row.sort()
		# tmpAr = [(index,rank)]
		tmpAr = [ (row[i][1],i) for i in range(len(row)) ]
		tmpAr.sort()
		for i in range(len(tmpAr)):
			outMatrix[i][j] = tmpAr[i][1]
	#print( outMatrix )
	return outMatrix

def readFlatFile( fileStr ):

	inFile = open( fileStr, 'r' )
	outAr = []
	for line in inFile:
		if line.startswith( '#' ):
			continue
		line = line.rstrip()
		outAr += [line]
	inFile.close()
	return outAr

def parseHeader( headerLine, samplesAr ):

	lineAr = headerLine.split( '\t' )
	# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
	# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
	# (9+4i) sample_i_FPKM (10+4i) sample_i_conf_lo (11+4i) sample_i_conf_hi
	# (12+4i) sample_i_status
	
	geneIndex = lineAr.index( 'gene_short_name' )
	if geneIndex == -1:
		print( 'ERROR: gene_short_name not found in FPKM file header' )
		exit()
	indexAr = [None] * len( samplesAr )
	for i in range( len(lineAr) ):
		if lineAr[i].endswith( 'FPKM' ):
			s1 = lineAr[i][:(lineAr[i].find('_'))]
			s2 = lineAr[i][:(lineAr[i].rfind('_'))]
			try:
				ind1 = samplesAr.index(s1)
			except ValueError:
				ind1 = -1
			try:
				ind2 = samplesAr.index(s2)
			except ValueError:
				ind2 = -1
			if ind1 != -1:
				indexAr[ind1] = i
			elif ind2 != -1:
				indexAr[ind2] = i
	if None in indexAr:
		print( 'ERROR: not all samples are in FPKM file' )
		exit()
	return indexAr, geneIndex

def writeOutput( outFileStr, sampleNamesAr, outMatrix, geneAr, baseJGenes, h3vGenes, info ):
	
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' )
	header = 'gene' 
	cC = 65
	for n in sampleNamesAr :	#+ ','.join( sampleNamesAr )
		header+= '\t{:c}\t{:s}\t{:c}'.format( cC, n, cC+1 )
		cC+= 2
	if baseJGenes != None:
		header += '\tbaseJ'
	if h3vGenes != None:
		header += '\th3v'
	outFile.write( header + '\n' )
	
	# loop through genes
	for i in range( len(geneAr) ):
		outStr = geneAr[i] + '\t'
		outStrAr = []
		# loop through samples
		for x in outMatrix[i]:
			outStrAr += ['{:.5f}'.format(x) ]*3
		outStr += '\t'.join( outStrAr )
		if baseJGenes != None:
			if geneAr[i] in baseJGenes:
				outStr += '\tNA'
			else:
				outStr += '\t0'
		if h3vGenes != None:
			if geneAr[i] in h3vGenes:
				outStr += '\tNA'
			else:
				outStr += '\t0'
		outFile.write( outStr + '\n' )
	outFile.close()

def processInputs( geneFileStr, fpkmFileStr, sampleNamesAr, baseJFileStr, h3vFileStr, absThreshold ):
	info = "#from_script:heatmap_rna_foldchange.py;"
	geneAr = readGeneFile( geneFileStr )
	print( 'Read gene file' )
	gind = geneFileStr.rfind('.')
	if absThreshold == None:
		outFileStr = geneFileStr[:gind] + '_heatmap_output_{:d}.tsv'.format( int( PERCENTILE*100 ))
		info += "percentile:{:d};".format( int(PERCENTILE*100))
	else:
		outFileStr = geneFileStr[:gind] + '_heatmap_output_{:d}.tsv'.format( int(absThreshold))
		info += "threshold:{:d};".format( int(absThreshold) )
	outMatrix = readFPKMFile( fpkmFileStr, sampleNamesAr, geneAr )
	print( 'Read FPKM file' )
	changeMatrix, numList = computeFoldChange( outMatrix )
	print( 'Calculated fold changes' )
	if absThreshold == None:
		thres = calculatePercentile( PERCENTILE, numList )
	else:
		thres = absThreshold
	correctMatrix = correctFoldChange( changeMatrix, thres )
	print( 'Corrected fold changes' )
	#orderedMatrix, outGenes = orderGenes( correctMatrix, geneAr )
	#print( 'Genes clustered' )
	#print( changeMatrix )
	info += 'gene_file:{:s};fpkm_file:{:s};'.format( os.path.basename( geneFileStr), os.path.basename( fpkmFileStr ) )
	if baseJFileStr != None:
		baseJGenes = readFlatFile( baseJFileStr )
		info += 'baseJ_file:{:s};'.format( os.path.basename( baseJFileStr ) )
		print( 'Read base J gene file' )
	else:
		baseJGenes = None
	if h3vFileStr != None:
		h3vGenes = readFlatFile( h3vFileStr )
		info += 'h3v_file:{:s};'.format( os.path.basename( h3vFileStr ) )
		print( 'Read H3V gene file' )
	else:
		h3vGenes = None
	info += 'base_sample:{:s};'.format( sampleNamesAr[0] )
	if len( sampleNamesAr ) == 2:
		info += 'compare_sample:{:s}'.format( sampleNamesAr[1] )
	elif len( sampleNamesAr ) > 2:
		info += 'compare_samples:{:s}'.format( ','.join( sampleNamesAr[1:] ) )
	writeOutput( outFileStr, sampleNamesAr[1:], correctMatrix, geneAr, baseJGenes, h3vGenes, info )
	print( 'Done' )
	
def parseInputs( argv ):
	baseJFile = None
	h3vFile = None
	absThreshold = None
	indStart = 0
	# get baseJ and/or H3V0 files and/or absolute threshold
	for i in range( 3 ):
		if argv[i].startswith( '-j=' ):
			baseJFile = argv[i].replace('-j=','')
			indStart += 1
		elif argv[i].startswith( '-h=' ):
			h3vFile = argv[i].replace('-h=','')
			indStart += 1
		elif argv[i].startswith( '-t=' ):
			try:
				absThreshold = float( argv[i].replace('-t=','') )
				indStart += 1
			except ValueError:
				print( 'ERROR: absolute threshold must be numeric' )
				exit()
	
	geneFileStr = argv[indStart]
	fpkmFileStr = argv[indStart+1]
	sampleNamesAr = []
	for i in range( indStart+2, len(argv) ):
		sampleNamesAr += [ argv[i] ]
	processInputs( geneFileStr, fpkmFileStr, sampleNamesAr, baseJFile, h3vFile, absThreshold )

if __name__ == "__main__":
	if len(sys.argv) < 5 :
		print ("Usage: python3.4 heatmap_rna_foldchange.py [-h=h3vo_gene_file] [-j=base_j_gene_file] [-t=threshold] <gene_file> <combined_fpkm_file> <base_name> <compare_name> [compare_name]*")
	else:
		parseInputs( sys.argv[1:] )
