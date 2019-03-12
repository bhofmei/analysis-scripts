import sys, math, glob, multiprocessing, subprocess, os, bisect

# Usage: python3.4 heatmap_from_bed_allc_pe.py [-1|-2] [-a|-r|-g|-n] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] [-f=ortholog_file] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95

ORTHO_FILE='/Users/bhofmeister/Documents/Research/collaborations/cmt3/arabidopsis_compare/Esalsugineum_Athaliana_bin.out'
#ORTHO_FILE='/escratch4/bth29393/bth29393_Jun_02/Ath_Esa_pairs.txt'

def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, genesIncluded, orthoFile, speciesInt ):
	''' 
	'''
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nGenes included: {:s}\nSpecies: {:d}\nPercentile for correction: {:.3f}\nSamples included: {:s}\n'.format( numBins, numBinsStream, streamSize, genesIncluded, speciesInt, percentile, ', '.join(sampleNamesAr) ) )
	if genesIncluded == 'all':
		subsetAr = None
	else:
		print( 'Reading A. thaliana - E. salsugineum ortholog file...' )
		subsetAr = readGBMPairs( genesIncluded, orthoFile, speciesInt )
	
	print( 'Reading FPKM file...' )
	fpkmAr =  getFPKM( fpkmFileStr, subsetAr )
	print( 'LEN:', len(fpkmAr) )
	print( 'Reading GFF file...' )
	geneDict, chrmFormat = readGFF( gffFileStr )
	if genesIncluded == 'all':
		outFileStrAr = [ '{:s}_{:d}_{:d}_{:d}.csv'.format( sampleName, streamSize, numBins, int(percentile*100) ) for sampleName in sampleNamesAr ]
	else:
		outFileStrAr = [ '{:s}_{:d}_{:d}_{:d}_{:s}.csv'.format( sampleName, streamSize, numBins, int(percentile*100), genesIncluded ) for sampleName in sampleNamesAr ]
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	results = [ pool.apply_async( processFile, args=(bedFileStrAr[i], outFileStrAr[i], fpkmAr, geneDict, numBins, numBinsStream, streamSize, percentile, chrmFormat )) for i in range( len(bedFileStrAr) ) ]
	fin = [ p.get() for p in results ]
	if sum(fin) != len( bedFileStrAr ):
		print( 'ERROR: not enough files written' )
		exit()
	print( 'Done' )

def readGBMPairs( genesIncluded, orthoFile, speciesInt ):
	
	outAr = []
	inFile = open( orthoFile, 'r' )
	ind = speciesInt * 2 + -2
	for line in inFile:
		lineAr = line.rstrip().split()
		i = lineAr[0].rfind('.')
		if i != -1:
			lineAr[0] = lineAr[0][:i]
		# (0) athaliana cds (1) UM/CG (2) eutrema cds (3) UM/CG
		try:
			# ortholog
			if genesIncluded == 'orthologs':
				bisect.insort( outAr, lineAr[ind] )
			# gbm
			elif genesIncluded == 'gBM_orthologs' and lineAr[1] == 'CG':
				bisect.insort( outAr, lineAr[ind] )
			# non-gBM
			elif genesIncluded == 'non_gBM_orthologs' and lineAr[1] == 'UM':
				bisect.insort( outAr, lineAr[ind] )
		except IndexError:
			print( 'ERROR:', line )
	inFile.close()
	return outAr

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

def getFPKM( fpkmFileStr, subsetAr ):
	'''
		returns the array or genes ordered by fpkm
	'''
	# check if exists
	checkFPKM( fpkmFileStr )
	
	fpkmFile = open( fpkmFileStr + ".txt", 'r' )
	fpkmAr = []
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		i = lineAr[0].rfind('.')
		if i != -1:
			lineAr[0] = lineAr[0][:i]
		if subsetAr == None or bisectIndex( subsetAr, lineAr[0] ) != None:
			fpkmAr += [lineAr[0]]
	fpkmFile.close()
	return fpkmAr

def checkFPKM( fpkmFileStr ):
	'''
		checks of the self-created tab-deliminated fpkm file exists
		it it exists but is outdated, it is recreated
		if it doesn't exist, it is created
	'''
	
	fileExists = os.path.isfile( fpkmFileStr + ".txt" )
	
	# if exists - check date modified and file size
	if fileExists:
		mTime = os.path.getmtime( fpkmFileStr )
		tTime = os.path.getmtime( fpkmFileStr + ".txt" )
		fSize = os.path.getsize( fpkmFileStr + ".txt" )
		if tTime < mTime or fSize < 100:
			print( 'Updating {:s}...'.format(fpkmFileStr + ".txt") )
			createFPKMTextFile( fpkmFileStr )
	# doesn't exist - create
	else:
		createFPKMTextFile( fpkmFileStr )
		print( 'Creating {:s}...'.format(fpkmFileStr + ".txt") )

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

def readFPKM( fpkmFileStr ):
	'''
		creates a dictionary for genes by fpkm value
		fpkm value is key which points to an array of gene names
		return fpkm dictionary
	'''
	fpkmFile = open( fpkmFileStr, 'r' )
	fpkmDict = {}
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		# Adam's output
		# (0) locus (1) coverage (2) FPKM
		#print( len( lineAr ) )
		if len( lineAr ) == 3:
			# header
			if line.startswith( 'locus' ):
				continue
			fpkm = float( lineAr[2] )
			name = lineAr[0]
			if fpkmDict.get( fpkm ) == None:
				fpkmDict[fpkm] = [name]
			# case 2: fpkm in dict
			else:
				fpkmDict[fpkm] += [name]
		# Cufflinks output
		# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
		# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
		# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
		else:
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
					# case 2: fpkm in dict
					else:
						fpkmDict[fpkm] += [n]
	return fpkmDict

def readGFF(gffFileStr):
	
	gffFile = open (gffFileStr, 'r' )
	gffDict = {}
	chrmFormat = None
	
	for line in gffFile:
		line = line[:-1]
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		chrm = lineAr[0]
		if chrmFormat == None:
			if chrm.isdigit():
				chrmFormat = 'digit'
			elif chrm.startswith( 'Chr0' ):
				chrmFormat = 'zeroed'
			elif chrm.startswith( 'Chr' ):
				chrmFormat = 'chr'
			else:
				chrmFormat = 'asis'
		strand = lineAr[6]
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			# put into geneDict
			gffDict[name] = (chrm, start, end, strand)
	
	gffFile.close()
	return gffDict, chrmFormat

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	n = notesStr[adIndex:]
	if endIndex != -1:
		n = notesStr[adIndex:endIndex+adIndex]
	
	pInd = n.rfind( '.' )
	if pInd != -1:
		return n[:pInd]
	else:
		return n

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		if fileStr.endswith( '.bed' ):
			leftIndex = fileStr.rfind('/')
			rightIndex = fileStr.rfind('.')
			sampleName = fileStr[leftIndex+1:rightIndex]
			sampleNamesAr += [ sampleName ]
		else:
			# allC file
			leftIndex = fileStr.rfind('/')
			sampleName = fileStr[leftIndex+1:]
			aIndex = sampleName.find( '_' )
			bIndex = sampleName[aIndex+1:].find( '_' )
			sampleNamesAr += [ sampleName[:bIndex+aIndex+1]+'_mCG' ]
	return sampleNamesAr

def processFile( fileStr, outFileStr, fpkmAr, gffDict, numBins, numBinsStream, streamSize, percentile, chrmFormat ):
	isBed = checkFile( fileStr )
	if isBed:
		print( 'Reading {:s}...'.format( fileStr ) )
		bedDict = readBed( fileStr )
		print( 'Processing {:s}'.format( fileStr ) )
		outMatrix, numList = processBed( bedDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize )
		thresh = calculatePercentile( percentile, numList )
		writeOutput( outFileStr, outMatrix, thresh )
		print( 'Output written for {:s}'.format( outFileStr ) )
		return 1
	else:
		print( 'Reading {:s}...'.format( fileStr ) )
		allCDict = readAllC( fileStr, chrmFormat )
		print( 'Processing {:s}'.format( fileStr ) )
		outMatrix, numList = processAllC( allCDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize )
		thresh = calculatePercentile( percentile, numList )
		writeOutput( outFileStr, outMatrix, thresh )
		print( 'Output written for {:s}'.format( outFileStr ) )
		return 1

def checkFile( fileStr ):
	if fileStr.endswith( '.bed' ):
		return True
	else:
		# allC file
		return False

def readBed( bedFileStr ):
	'''
		creates a dictionary with each scaffold and those dictionaries are a 
		dictionary
		for the positions of coordinates with frequency count from the bed file
		{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
		return the dictionary
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		lineAr =line.rstrip().split('\t')
		try:
			curScaf = lineAr[0]
			pos = int( lineAr[1] ) + 1
		except ValueError:
			pass
		else:
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
	return bedDict

def processBed( bedDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize ):
	# fpkmDict organized {fpkm: (gene1, gene2, ...)}
	# geneDict organized {name: (scaffold, start, end, strand)}
	countGenes = 0
	outMatrix = []
	numList = []
	# loop by gene length
	for gene in fpkmAr:
		info = geneDict.get( gene )
		# not actually a gene - don't count it
		if info == None:
			continue
		countGenes += 1
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = varByRegion(bedDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = varByRegion(bedDict, chrm, end, downstream, numBinsStream)
		# gene
		curGeneVarAr = varByRegion(bedDict, chrm, start, end, numBins)
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				outAr = curUpVarAr + curGeneVarAr + curDownVarAr
			else:
				outAr = curGeneVarAr
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				outAr = curDownVarAr + curGeneVarAr + curUpVarAr
			else:
				curGeneVarAr.reverse()
				outAr = curGeneVarAr
		# add to numList
		numList += outAr
		outMatrix.append( outAr )
	print( 'num genes included:', countGenes )
	return outMatrix, numList

def varByRegion( bedDict, chrm, start, end, numBins ):
	''' 
		takes in the variant dictionary generated for a gene and the start and 
		end positions of that gene. Note the dictionary is set up so the key is 	
		the position
		on the scaffold and the value is the frequency of that variant
		returns one array with number of variants separated by bin
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print 'binWidth', binWidth
	varAr = [0] * numBins
	
	# loop for each bin
	for bin in range(numBins - 1):
		# loop for each position in that bin
		for pos in range(binWidth):
			key = bin * binWidth + pos + start
			try:
				dictEntry = bedDict.get(chrm).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					varAr[bin] += dictEntry
			except AttributeError:
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = bedDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				varAr[-1] += dictEntry
		except AttributeError:
			pass
	
	varAr = adjustArrayLength( binWidth, varAr )
	return varAr

def adjustArrayLength( length, inAr ):
	outAr = [0] * len( inAr )
	for i in range( len(inAr) ):
		outAr[i] = float( inAr[i] ) / float( length )
	return outAr

def calculatePercentile( percentile, numList ):
	
	if percentile == 1:
		return max( numList )
	numList.sort()
	ind = math.ceil( percentile * len( numList ) - 1 )
	try:
		p = numList[ind]
	except IndexError:
		return numList[-1]
	if p == 0:
		print( '***** not percentile corrected *****' )
		return numList[-1]
	return p

def readAllC( allCFileStr, chrmFormat ):
	'''
		This allC file should be the allC information for all chromosomes in one 
		file not just a single chromosome
	'''
	allCFile = open( allCFileStr, 'r' )
	allCDict = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if lineAr[3].startswith( 'CG' ):
			# no dictionary for scaffold
			chrm = lineAr[0]
			# handle formatting
			if chrmFormat == 'digit' and chrm.isdigit() == False:
				chrm = chrm.replace( 'Chr', '' )
			elif chrmFormat == 'zeroed' and chrm.isdigit():
				chrm = 'Chr{:02d}'.format( int(chrm) )
			elif chrmFormat == 'chr' and chrm.isdigit():
				chrm = 'Chr{:d}'.format( int( chrm ) )
			if allCDict.get( chrm ) == None:
				allCDict[ chrm ] = {}
			allCDict[chrm][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
	
	allCFile.close()
	return allCDict

def processAllC( allCDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize ):
	
	outMatrix = []
	numList = []
	#print( allCDict.keys() )
	# loop by gene length
	for gene in fpkmAr:
		info = geneDict.get( gene )
		if info == None:
			continue
		#print( info )
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = methByRegion(allCDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = methByRegion(allCDict, chrm, end, downstream, numBinsStream)
		
		curGeneVarAr = methByRegion(allCDict, chrm, start, end, numBins)
		#print( 'curRepeatVarAr', curRepeatVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				outAr = curUpVarAr + curGeneVarAr + curDownVarAr
			else:
				outAr = curGeneVarAr
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				outAr = curDownVarAr + curGeneVarAr + curUpVarAr
			else:
				curGeneVarAr.reverse()
				outAr = curGeneVarAr
		outMatrix.append( outAr )
		numList += outAr
	return outMatrix, numList

def methByRegion(allCDict, chrm, start, end, numBins):
	''' 
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
	#print ('binWidth', binWidth)
	methAr = [0] * numBins
	totalAr = [0] * numBins
	
	# loop for each bin
	for bin in range(numBins - 1):
		# loop for each position in that bin
		for pos in range(binWidth):
			key = bin * binWidth + pos + start
			#print( key )
			try:
				dictEntry = allCDict.get(chrm).get(key)
				# only want to include sites we have info for
				if dictEntry != None:
					methAr[bin]+= dictEntry[0]
					totalAr[bin] += dictEntry[1]
			# chrm not there
			except AttributeError:
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = allCDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				methAr[-1] += dictEntry[0]
				totalAr[-1] += dictEntry[1]
		except AttributeError:
			pass
	z = zip( methAr, totalAr )
	methAr = [ (m,t) for m,t in z ]
	outAr = calculateMethylation( methAr, binWidth )
	return outAr

def calculateMethylation( methAr, length ):
	
	outAr = [-1] * len( methAr )
	for i in range(len(methAr)) :
		if methAr[i][0] == 0 and methAr[i][1] == 0:
			outAr[i] = 0
		elif methAr[i][0] != 0 and methAr[i][1] == 0:
			print( methAr[i] )
		else:
			outAr[i] = (float(methAr[i][0]) / float(methAr[i][1])) / length
	return outAr

def writeOutput( outFileStr, outMatrix, threshold ):

	outFile = open( outFileStr, 'w' )
	geneCount = len( outMatrix )
	
	for gene in outMatrix:
		for i in range( 0, len(gene) ):
			if gene[i] > threshold:
				gene[i] = threshold
		outStr = 'G{:08},'.format( geneCount )
		outStrAr = [ '{:.5f}'.format(x) for x in gene ]
		outStr += ','.join( outStrAr )
		outFile.write( outStr + '\n' )
		geneCount -= 1
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	percentile = PERCENTILE
	genesIncluded = 'all'
	speciesInt = 0
	startInd = 0
	orthoFile = ORTHO_FILE
	
	for i in range( min(11, len(argv)-3 ) ):
		# number of processors
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-s=' ):
			try:
				streamSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: upstream/downstream size must be integer' )
				exit()
		elif argv[i].startswith( '-b=' ):
			try:
				numBinsAr = argv[i][3:].split(',')
				numBins = int( numBinsAr[0] )
				if len( numBinsAr ) == 2:
					numBinsStream = int( numBinsAr[1] )
				else:
					numBinsStream = int( numBinsAr[0] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of bins must be integer' )
				exit()
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		elif argv[i] == '-r':
			# check -g or -n not set
			if genesIncluded in ['non_gBM_orthologs','gBM_orthologs']:
				print('ERROR: can only specify one of -r, -g, -n options' )
				exit()
			genesIncluded = 'orthologs'
			startInd += 1
		elif argv[i] == '-g':
			# check -n or -r not set
			if genesIncluded in ['non_gBM_orthologs','orthologs']:
				print('ERROR: can only specify one of -r, -g, -n options' )
				exit()
			genesIncluded = 'gBM_orthologs'
			startInd += 1
		elif argv[i] == '-n':
			# check -g or -r not set
			if genesIncluded in ['gBM_orthologs','orthologs']:
				print('ERROR: can only specify one of -r, -g, -n options' )
				exit()
			genesIncluded = 'non_gBM_orthologs'
			startInd += 1
		elif argv[i] == '-a':
			if genesIncluded in ['gBM_orthologs','non_gBM_orthologs','orthologs']:
				print( 'ERROR: cannot specify -a with -r, -n, or -g' )
				exit()
			genesIncluded = 'all'
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			orthoFile = argv[i][3:]
			startInd += 1
		elif argv[i] == '-1':
			if speciesInt != 0:
				print( 'ERROR: cannot specify -1 with -2' )
				exit
			speciesInt = 1
			startInd += 1
		elif argv[i] == '-2':
			if speciesInt != 0:
				print( 'ERROR: cannot specify -2 with -1' )
				exit()
			speciesInt = 2
			startInd += 1
	# end for
	
	if genesIncluded == 'all' and speciesInt != 0:
		print( 'ERROR: do not specify species int when using all genes' )
		exit()
	elif speciesInt == 0 and genesInclude != 'all':
		speciesInt = 2
	
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, genesIncluded, orthoFile, speciesInt )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 heatmap_from_bed_allc_pe.py [-1|-2] [-a|-r|-g|-n] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] [-f=ortholog_file] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*\n-1\tuse genes listed first in ortholog file\n-2\tuse genes listed second in ortholog file\n-a\tinclude all genes\n-r\tinclude only orthologs\n-g\tinclude only gBM orthologs\n-n\tinclude only non-gBM orthologs\n-f=ortholog_file\tuse this ortholog file instead\np=num_proc\tnumber of processors to use [default 2]\n-b=N,W\tuse N bings for gene body and W bins for upstream/downstream\n-s\tnumber of base pairs to include upstream and downstream\n-t=D\tpercentile to correct largest values by")
	else:
		parseInputs( sys.argv[1:] )
