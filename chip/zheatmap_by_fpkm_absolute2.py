import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 heatmap_by_fpkm_absolute.py <gff_file> <fpkm_file> <distance> <bed_file> [bed_files]* [percentile_correction]
BINWIDTH=50

def readGFF(gffFileStr):
	
	gffFile = open (gffFileStr, 'r' )
	geneDict = {}
	
	for line in gffFile:
		line = line[:-1]
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		scaffold = lineAr[0]
		strand = lineAr[6]
		#print ( "GFF:",lineAr)
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			# put into geneDict
			geneDict[name] = (scaffold, start, end, strand)
	
	gffFile.close()
	#print( geneDict )
	return geneDict

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]
	
def checkFPKM( fpkmFileStr ):
	'''
		checks of the self-created tab-deliminated fpkm file exists
		it it exists but is outdated, it is recreated
		if it doesn't exist, it is created
	'''
	
	fileExists = os.path.isfile( fpkmFileStr + ".txt" )
	
	# if exists - check date modified
	if fileExists:
		mTime = os.path.getmtime( fpkmFileStr )
		tTime = os.path.getmtime( fpkmFileStr + ".txt" )
		if tTime < mTime:
			print( 'Updating {:s}...'.format(fpkmFileStr + ".txt") )
			createFPKMTextFile( fpkmFileStr )
	# doesn't exist - create
	else:
		createFPKMTextFile( fpkmFileStr )
		print( 'Creating {:s}...'.format(fpkmFileStr + ".txt") )
	
def readFPKM( fpkmFileStr ):
	'''
		creates a dictionary for genes by fpkm value
		fpkm value is key which points to an array of gene names
		return fpkm dictionary
	'''
	fpkmFile = open( fpkmFileStr, 'r' )
	
	# headers for fpkm
	# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id (4) gene_short_name
	# (5) tss_id (6) locus (7) length (8) coverage (9) FPKM (10) FPKM_conf_low
	# (11) FPKM_conf_high (12) FPKM_status
	fpkmDict = {}
	
	for line in fpkmFile:
		line = line[:-1]
		lineAr = line.split('\t')
		if lineAr[12] == 'OK':
			fpkm = float( lineAr[9] )
			name = lineAr[4]
			if name == '-':
				continue
			nameAr = name.split(',')
			for n in nameAr:
				# case 1: fpkm not in fpkmDict
				if fpkmDict.get( fpkm ) == None:
					fpkmDict[fpkm] = [n]
				# case 2: fpkm in fpkmDict
				else:
					fpkmDict[fpkm] += [n]
	return fpkmDict

def createFPKMTextFile( fpkmFileStr ):
	'''
		creates the tab-deliminated FPKM file
		useful so gene order is the same across samples
	'''
	# fpkmDict set up as {fpkm:[gene1,gene2,...], fpkm2:[gene4],...}
	fpkmOutFile = open( fpkmFileStr + ".txt", 'w' )
	fpkmDict = readFPKM( fpkmFileStr )
	
	for key in sorted(fpkmDict.keys(), reverse=True):
		genes = fpkmDict[key]
		# for each gene
		for gene in genes:
			fpkmOutFile.write( "{:s}\t{:.4f}\n".format( gene, key ) )
			
	fpkmOutFile.close()

def getFPKM( fpkmFileStr ):
	'''
		returns the array of genes ordered by fpkm and count of genes
	'''
	# check if exists
	checkFPKM( fpkmFileStr )
	
	fpkmFile = open( fpkmFileStr + ".txt", 'r' )
	fpkmGeneAr = []
	#fpkmAr = []
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		fpkmGeneAr += [lineAr[0]]
	fpkmFile.close()
	return fpkmGeneAr

def fullProcessBed( bedFileStr, gffDict, fpkmAr, distance, percentile ):
	
	print( 'starting to process', bedFileStr )
	sampleName = getSampleName( bedFileStr )
	bedDict = readBed( bedFileStr )
	print( 'bed read for', sampleName )
	outMatrix, numList, geneList = processBed( gffDict, distance, bedDict, fpkmAr )
	#print( outMatrix )
	threshold = calculatePercentile( percentile, numList )
	if percentile != 1:
		outFileStr = '{:s}_{:d}_{:d}.csv'.format( sampleName, distance, int( math.ceil(percentile*100)) )
	else:
		outFileStr = '{:s}_{:d}.csv'.format( sampleName, distance )
	writeOutput( outFileStr, outMatrix, geneList, threshold )
	print( 'output written for', sampleName )

def readBed(bedFileStr):
	'''
		creates a dictionary with each scaffold and those dictionaries are a 
		dictionary for the positions of coordinates with frequency count from 
		the bed file
		{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
		return the dictionary
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		line = line.replace('\n','')
		lineAr =line.split('\t')
		try:
			curScaf = lineAr[0]
			pos = int( lineAr[1] )
		except ValueError:
			pass
		else:
			#print lineAr
			# no dictionary for scaffold
			if bedDict.get(curScaf) == None:
				bedDict[curScaf] = {pos:1}
			# dictionary for scaffold but position not included
			elif bedDict.get(curScaf).get(pos) == None:
				bedDict[curScaf][pos] = 1
			# dictionary for scaffold and position included
			else:
				bedDict[curScaf][pos] += 1
	bedFile.close()
	#print( "--bed read -- dictionary created" )
	return bedDict

def processBed( gffDict, distance, bedDict, fpkmAr ):
	# geneDict organized {name: (scaffold, start, end, strand)}
	# bedDict organized {scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
	
	numList = []
	outMatrix = []
	geneList = []
	
	for i in range( len(fpkmAr) ):
		info = gffDict.get( fpkmAr[i] )
		if info == None:
			continue
		scaffold = info[0]
		bedDictS = bedDict.get( scaffold )
		if bedDictS == None:
			continue
		start = info[1]
		end = info[2]
		strand = info[3]
		strandNum = (0 if strand == '+' else 1 )
		curArF, curArT = bedCountRegion( bedDictS, start, end, distance, strandNum )
		numList += [ x for x in curArF if x != -1 ] + [ x for x in curArT if x != -1 ]
		outMatrix.append( curArF + [0]*2 + curArT )
		geneList += fpkmAr[i]
	return outMatrix, numList, geneList

def bedCountRegion( bedDict, start, end, distance, strandNum ):
	numBins = int(distance//BINWIDTH)
	outArF = [-1] * (numBins*2)
	outArT = [-1] * (numBins*2)
	geneLength = end - start + 1
	
	# long enough genes
	if geneLength >= 2*distance:
		for bin in range( len( outArF ) ):
			for i in range( BINWIDTH ):
				# 5` for + genes
				if outArF[bin] == -1:
					outArF[bin] = 0
				dictEntry = bedDict.get((start - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArF[bin] += dictEntry
				# 3` for + genes
				if outArT[bin] == -1:
					outArT[bin] = 0
				dictEntry = bedDict.get((end - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArT[bin] += dictEntry
	else:
		indStop = int((geneLength / 2 ) // BINWIDTH)
		ind1 = 0
		ind2 = numBins  + indStop
		ind3 = numBins - indStop
		ind4 = len( outArT )
		
			
		# 5` for + genes (3` for - genes)
		for bin in range( ind1, ind2 ):
			for i in range( BINWIDTH ):
				if outArF[bin] == -1:
					outArF[bin] = 0
				dictEntry = bedDict.get((start - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArF[bin] += dictEntry
		# 3` for + genes (5` for - genes)
		for bin in range( ind3, ind4 ):
			for i in range( BINWIDTH ):
				if outArT[bin] == -1:
					outArT[bin] = 0
				dictEntry = bedDict.get((end - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArT[bin] += dictEntry
	# normalize by distance
	outArF = [ float(x)/BINWIDTH for x in outArF if x != -1]
	outArT = [ float(x)/BINWIDTH for x in outArT if x != -1 ]
	# need to reverse for negative strand and return opposite
	if strandNum == 1:
		outArF.reverse()
		outArT.reverse()
		return outArT, outArF
	else:
		return outArF, outArT

	
def writeOutput( outFileStr, outMatrix, geneList, threshold ):
	
	outFile = open( outFileStr, 'w' )
	# loop through gene
	for i in range( len(outMatrix) ):
		geneID = 'G{:08d}'.format(len( outMatrix ) - i)
		if max( outMatrix[i] ) > threshold:
			outMatrix[i] = [ ( threshold if x > threshold else x ) for x in outMatrix[i] ]
		arStr = [ ('NA' if x == -1 else '{:.5f}'.format(x)) for x in outMatrix[i] ]
		outStr = geneID + ',' + geneList[i] + ',' + ','.join( arStr ) + '\n'
		outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	try:
		percentile = float( argv[0] )
		if percentile > 1:
			percentile = percentile / 100
		startL = 1
	except ValueError:
		percentile = 1
		startL = 0
	gffFileStr = argv[startL]
	fpkmFileStr = argv[startL+1]
	try:
		distance = int( argv[startL+2] )
		if distance % 100 != 0:
			print( 'ERROR: distance must be divisible by 100' )
			exit()
	except ValueError:
		print( 'ERROR: distance must be an integer' )
		exit()
	bedFileStrAr =[]
	for i in range( startL+3, len(argv) ):
		bedFileStrAr += [ argv[i] ]
	processInputs( gffFileStr, fpkmFileStr, distance, percentile, bedFileStrAr )
		
def processInputs( gffFileStr, fpkmFileStr, distance, percentile, bedFileStrAr ):
	gffDict = readGFF( gffFileStr )
	print( 'read gff' )
	fpkmAr = getFPKM( fpkmFileStr )
	print( 'read fpkm' )
	
	#pool = multiprocessing.Pool( processes=3 )
	#t = [pool.apply( fullProcessBed, args=(f, gffDict, fpkmAr, distance ) ) for f in bedFileStrAr ]
	for f in bedFileStrAr:
		fullProcessBed( f, gffDict, fpkmAr, distance, percentile )
	print( 'Done' )

def getSampleName( fileStr):
	lind = fileStr.rfind('/')
	rind = fileStr.rfind('.')
	name = fileStr[lind+1:rind]
	return name

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
	

if __name__ == "__main__":
	if len(sys.argv) < 5 :
		print ("Usage: python3.4 heatmap_by_fpkm_absolute.py [percentile_correction] <gff_file> <fpkm_file> <distance> <bed_file> [bed_files]*")
	else:
		parseInputs( sys.argv[1:] )
