import sys, math, glob, multiprocessing, subprocess, os, random

# Usage: python3.4 chip_ortholog_compare2.py [-o=outID] [-p=num_proc] [-r=intergenic_region_width] [-n=num_features] [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id] [-i1=intergenic_file_species1] [-i2=<intergenic_species2] -a1=<allc_species1> -a2=<allc_species2> -b1<chip_bed_species1> -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>
# Using methylation information and gene orthologs creates heatmaps showing ChIP enrichment for two samples; includes intergenic and genic regions

# allC_files chip_bed_files ortholog_file gff_files intergenic_bed_files
# fasta_files

REGIONWIDTH=1000
NUMSAMPLES=1000
NUMCGS=40
NUMBINS=20
NUMBINSSTREAM=0
STREAMSIZE=0
PERCENTILE=0.95
NUMPROC=1
NUMINTERGENICRAND=25000

# A. intergenic regions
# make new list of regions
# break into 1 kb sections
# weighted methylation of those to rank
# 1. intergenic unmethylated regions
# 2. intergenic methylated regions (repeats/transposons?)
# B. genes
# 3. arabidopsis unmethylated genes and orthologs
# 4. arabidopsis methylated genes and orthologs
# 5. output

def processInputs( allcFileAr, chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr, intergenicBedAr, speciesIDAr, outID, numProc, regionWidth, numFeatures, numCGs, numBins, numBinsStream, streamSize, percentile ):
	# fpkmFileStr, orthologFileStr, numFeatures, outID, speciesIDAr, numProc
	# gffFileAr, intergenicBedAr, allcFileAr, chipBedAr, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile
	
	gffStr = [ os.path.basename( x ) for x in gffFileAr ]
	allcStr = [ os.path.basename( x ) for x in allcFileAr ]
	chipBedStr = [ os.path.basename( x ) for x in chipBedAr ]
	# write meta information
	print( numFeatures, regionWidth, numCGs, numBins, numBinsStream, streamSize, os.path.basename(orthoFileStr), gffStr, allcStr, chipBedStr, os.path.basename(fpkmFileStr) )
	info = '#from_script:chip_ortholog_compare.py numFeatures:{:d}; regionWidth:{:d}; numCGs:{:d}; numBins:{:d}; numBinsStream:{:d}; streamSize:{:d}; ortho_file:{:s}; gff_files:{:s}; allc_files:{:s}; chip_bed_files:{:s}; fpkm_file:{:s}'.format( numFeatures, regionWidth, numCGs, numBins, numBinsStream, streamSize, os.path.basename(orthoFileStr), ','.join(gffStr), ','.join(allcStr), ','.join(chipBedStr), os.path.basename(fpkmFileStr) )
	
	# get fpkm
	fpkmAr, fpkmValAr = getFPKM( fpkmFileStr )
	
	# get orthologs
	print( 'Reading {:s}...'.format( orthoFileStr ) )
	methAr, unmethAr = readOrtholog( orthoFileStr )
	# get gene lists (random/methylation rank then order fpkm)
	print( 'Determining orthologous genic regions...' )
	methArs, fpkmMethAr = trimGenic( methAr, fpkmAr, fpkmValAr )
	unmethArs, fpkmUnmethAr = trimGenic( unmethAr, fpkmAr, fpkmValAr )
	
	outFileAr = [ '{:s}_{:s}.tsv'.format( outID, speciesID ) for speciesID in speciesIDAr ]
	
	# per species -> can be done in parallel
	# gffFileStr, intergenicBedStr, allcFileStr, chipBedStr, outFileStr, unmethAr, methAr, speciesID, regionWidth, numCGs, numFeatures, numBins, numBinsStream, streamSize, percentile
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing species with {:d} processes'.format( numProc ) )
	results = [ pool.apply_async( processSpecies, args=( gffFileAr[i], intergenicBedAr[i], allcFileAr[i], chipBedAr[i], fastaFileAr[i], outFileAr[i], unmethArs[i], methArs[i], speciesIDAr[i], regionWidth, numCGs, numFeatures, numBins, numBinsStream, streamSize, percentile, fpkmMethAr, fpkmUnmethAr, info ) ) for i in range(len(speciesIDAr)) ]
	metaMatrix = [ p.get() for p in results ]
	outMetaStr = '{:s}_meta.tsv'.format( outID )
	print( 'Writing metaplot output to {:s}...'.format( outMetaStr ) )
	writeMetaOutput( outMetaStr, metaMatrix, speciesIDAr, numBins, numBinsStream, info )
	print( 'Done.' )

def processSpecies( gffFileStr, intergenicBedStr, allcFileStr, chipBedStr, fastaFileStr, outFileStr, unmethAr, methAr, speciesID, regionWidth, numCGs, numFeatures, numBins, numBinsStream, streamSize, percentile, fpkmMethAr, fpkmUnmethAr, info ):

	# check intergenic
	if intergenicBedStr == None:
		ofStr = 'intergenic_{:s}.tsv'.format( speciesID )
		print( 'Creating {:s}...'.format( ofStr ) )
		readGFFIntergenic( gffFileStr, ofStr )
		intergenicBedStr = ofStr
	print( 'Preparing intergenic regions for species {:s}...'.format( speciesID ) )
	print( 'Reading {:s}...'.format( fastaFileStr ) )
	fastaDict = readFasta( fastaFileStr )
	print( 'Reading {:s}...'.format( intergenicBedStr ) )
	intergenicAr = readIntergenic( intergenicBedStr, regionWidth )
	print( 'Trimming...' )
	trimIntergenicAr = trimIntergenic( intergenicAr, fastaDict, numCGs )
	print( 'Intergenic regions of {:s} trimmed from {:d} to {:d}'.format( speciesID, len(intergenicAr), len(trimIntergenicAr) ) )
	print( 'Reading {:s}...'.format( allcFileStr ) )
	allCDict = readAllC( allcFileStr )
	print( 'Calculating and ranking intergenic regions for {:s}...'.format( speciesID) )
	intergenicDict = calculateIntergenicRegions( trimIntergenicAr, allCDict )
	writeIntergenicValues( 'intergenic_values_{:s}.tsv'.format( speciesID ), intergenicDict )
	highestIntergenic, lowestIntergenic = rankIntergenicRegions( intergenicDict, numFeatures )
	
	# gff dict
	print( 'Reading {:s}...'.format( gffFileStr ) )
	gffDict = readGFF( gffFileStr )
	
	# process regions
	print( 'Processing all regions for {:s}...'.format( speciesID ) )
	if numBinsStream != 0 and streamSize == 0:
		numBinsStream = 0
		
	nameMatrix, outMatrix, numList, metaAr = processBedRegions( chipBedStr, lowestIntergenic, highestIntergenic, unmethAr, methAr, gffDict, numBins, numBinsStream, streamSize )
	# threshold
	thresh = calculatePercentile( percentile, numList )
	fpkmAr = ['NA']*len(lowestIntergenic) + ['NA'] * len(highestIntergenic) + fpkmMethAr + fpkmUnmethAr
	
	# write output
	print( 'Writing species {:s} output to {:s}'.format( speciesID, outFileStr ) )
	#info = '#numFeatures:{:d}_regionWidth:{:d}_numCGs:{:d}_numBins:{:d}_numBinsStream:{:d}_streamSize:{:d}_percentile:{:.3f}'.format( numFeatures, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile )
	
	writeOutput( outFileStr, nameMatrix, outMatrix, thresh, numBins, numBinsStream, info, fpkmAr)
	return metaAr

def readGFFIntergenic( gffFileStr, outFileStr ):
	'''
		reads the gene GFF to parse out intergenic regions
		intergenic is +/- 100bp between genes
		writes intergenic regions to file
	'''
	
	prevPos = 1
	curChrm = None
	
	gffFile = open( gffFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		
		# ignore non-genes
		if lineAr[2] != 'gene':
			continue
		
		chrm = lineAr[0]
		adStart = int( lineAr[3] ) - 100
		adEnd = int( lineAr[4] ) + 100
		
		# check new chrm
		if chrm != curChrm:
			#outStr += '{:s}\t{:d}\t{:d}\n'.format( chrm, 1, adStart )
			prevPos = adEnd
			curChrm = chrm
		# same chrm
		else:
			outFile.write('{:s}\t{:d}\t{:d}\n'.format( chrm, min(prevPos, adStart), max(prevPos, adStart) ) )
			prevPos = adEnd
	gffFile.close()

def readFasta( fastaFileStr ):
	''' 
		we will read in the whole fasta because it shouldn't be too big
		read into dictionary where chrm is the key and and the value is the
		sequence
		NOTE: we are making the sequence 1-based indexed and automatically 
		capitalizing all of the bases
		returns the created dictionary
	'''
	fastaFile = open( fastaFileStr, 'r' )
	fastaDict = {}
	chrm = None
	seq = ' '
	for line in fastaFile:
		line = line.rstrip()
		# chromosome headers
		if line.startswith('>'):
			# need to write old sequence
			if seq != ' ' and chrm != None:
				fastaDict[chrm] = seq
				seq = ' '
			lineAr = line.split(' ')
			chrm = lineAr[0].replace('>', '')
		
		# sequence
		else:
			seqL = line.upper()
			seq += seqL
			
	# handle last chrm read
	if chrm not in fastaDict:
		fastaDict[chrm] = seq
	fastaFile.close()
	return fastaDict

def readIntergenic( intergenicBedStr, regionWidth ):
	'''
		read intergenic file; as it reads a region, divides that into
		regionWidth sized subregions 
		return array of accpetable intergenic regions
		array contains tuples (chrm, subregion_start, subregion_end)
	'''
	
	intergenicAr = []
	
	intergenicBedFile = open( intergenicBedStr, 'r' )
	for line in intergenicBedFile:
		lineAr = line.rstrip().split( '\t' )
		chrm = lineAr[0]
		start = int(lineAr[1])
		end = int(lineAr[2])
		length = end - start + 1
		for i in range( length // regionWidth ):
			adStart = start + i*regionWidth
			adEnd = start + (i+1)*regionWidth - 1
			intergenicAr += [(chrm, adStart, adEnd )]
	intergenicBedFile.close()
	return intergenicAr

def trimIntergenic( intergenicAr, fastaDict, numCGs ):
	'''
		trim intergenic regions based on the number of CGs
		return array of intergenic regions with enough CGs
		array contains tuples (chrm, start, end, numCGs)
	'''
	
	outAr = []
	
	for region in intergenicAr:
		chrm = region[0]
		start = region[1]
		end = region[2]
		cCGs = fastaDict[chrm][start:end+1].count('CG')
		if cCGs >= numCGs:
			outAr += [ (chrm, start, end, cCGs) ]
	randAr = randomIntergenic( outAr )
	return randAr

def randomIntergenic( intergenicAr ):
	numFeatures = NUMINTERGENICRAND
	numRegions = len(intergenicAr)-1
	randIntSet = set()
	tmpAr = []
	# pick indexes
	while len( randIntSet ) < numFeatures:
		p = random.randint( 0, numRegions )
		if p in randIntSet:
			continue
		tmpAr += [ intergenicAr[p] ]
		randIntSet.add(p)
	return tmpAr

def readAllC( allCFileStr ):
	'''
		This allC file should be the allC information for all chromosomes in one 
		file not just a single chromosome
	'''
	mTypes = ['CG','CHG','CHH']
	allCFile = open( allCFileStr, 'r' )
	allCDict = {}
	for m in mTypes:
		allCDict[m] = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		mType = decodeMethType( lineAr[3] )
		if mType != False:
			# no dictionary for scaffold
			chrm = lineAr[0]
			if allCDict.get( chrm ) == None:
				allCDict[ chrm ] = {}
			allCDict[chrm][int(lineAr[1])] = ( mType, int(lineAr[4]), int( lineAr[5]) )
	
	allCFile.close()
	return allCDict

def decodeMethType( mStr ):
	
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return False
	else:
		return 'CHH'

def calculateIntergenicRegions( intergenicAr, allCDict ):
	'''
		goes through intergenic regions and computes weighted methylation
		puts regions in a dictionary based on weighted CG methylation
		dict key: weighted methylation, dict value: tuple (chrm, start, end, numCGs, weighted CHG, weighted CHH )
	'''
	intergenicDict = {}
	count = 0
	for region in intergenicAr:
		if count % 5000 == 0:
			print( 'intergenic region', count )
		count += 1
		# region: (chrm, start, end, numCGs )
		wMethCG, wMethCHG, wMethCHH = calculateWeightedMethylation( allCDict, region[0], region[1], region[2] )
		if wMethCG == -1:
			continue
		if intergenicDict.get( wMethCG ) == None:
			intergenicDict[wMethCG] = []
		intergenicDict[wMethCG] += [ (region[0], region[1], region[2], region[3], wMethCHG, wMethCHH) ]
	return intergenicDict
	
def calculateWeightedMethylation( allCDict, chrm, start, end ):
	'''
		compute weighted methylation of region in all 3 sequence contexts
		returns weighted methylation of CG, CHG, CHH; -1 if no reads for region
	'''
	readsMeth = [0,0,0]
	readsTotal = [0,0,0]
	mTypes = ['CG','CHG','CHH']
	#print( chrm, start, end )
	for pos in range( start, end+1 ):
		try:
			dictEntry = allCDict.get(chrm).get( pos )
			if dictEntry != None:
				mInd = mTypes.index( dictEntry[0] )
				readsMeth[mInd] += dictEntry[1]
				readsTotal[mInd] += dictEntry[2]
		except AttributeError:
			pass
	outAr = [ (-1 if readsMeth[i] == 0 and readsTotal[i] == 0 else float(readsMeth[i]) / float(readsTotal[i] ) ) for i in range(len(mTypes)) ]
	return outAr[0], outAr[1], outAr[2]

def rankIntergenicRegions( intergenicDict, numFeatures ):
	'''
		return list of highest methylated regions and lowest methylated regions
		for regions with same methylation, randomly orders/chooses regions
	'''
	highestAr = []
	lowestAr = []
	
	# get highest
	for key in sorted( intergenicDict.keys(), reverse=True ):
		regionAr = intergenicDict[key]
		shuf = list(range(len(regionAr)))
		random.shuffle( shuf )
		for i in shuf:
			highestAr += [ regionAr[i] ]
		if len( highestAr ) >= numFeatures:
			break
	# get lowest
	for key in sorted( intergenicDict.keys() ):
		regionAr = intergenicDict[key]
		shuf = list(range(len(regionAr)))
		random.shuffle( shuf )
		for i in shuf:
			lowestAr += [ regionAr[i] ]
		if len( lowestAr ) >= numFeatures:
			break
	# highest is ordered highest to lowest
	# lowest is ordered lowest to highest
	# make all of them ordered highest to lowest
	lowestAr = lowestAr[:numFeatures]
	lowestAr.reverse()
	highestAr = highestAr[:numFeatures]
	return highestAr, lowestAr

def writeIntergenicValues( outFileStr, intergenicDict ):
	
	
	outFile = open( outFileStr, 'w' )
	# header: chrm start end numCGs weightedMethylation
	outFile.write( 'chrm\tstart\tend\tweightedCGMethylation\tnumCGs\tweightedCHG\tweightedCHH\n' )
	for key in sorted( intergenicDict.keys() ):
		regionAr = intergenicDict[key]
		for region in regionAr:
			outFile.write( '{:s}\t{:d}\t{:d}\t{:.2f}\t{:d}\t{:.2f}\t{:.2f}\n'.format( region[0], region[1], region[2], key*100, region[3], region[4]*100, region[5]*100 ) )
	outFile.close()
	
def readGFF( gffFileStr ):
	
	gffFile = open (gffFileStr, 'r' )
	gffDict = {}
	
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
		strand = lineAr[6]
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			# put into geneDict
			gffDict[name] = (chrm, start, end, strand)
	
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
	
def readOrtholog( orthologFileStr ):
	
	methAr = []
	unmethAr = []
	isFirst = True
	
	orthoFile = open( orthologFileStr, 'r' )
	for line in orthoFile:
		if isFirst:
			isFirst = False
			continue
		lineAr = line.rstrip().split( )
		# (0) Athaliana_gene (1) gBM/UM (2) Esalsugineum_gene (3) gBM/UM
		atGene = lineAr[0][:-2]
		esGene = lineAr[2] + '.g'
		if lineAr[1] == 'CG':
			methAr += [ (atGene, esGene) ]
		elif lineAr[1] == 'UM' and lineAr[1] == 'UM':
			unmethAr += [ (atGene, esGene) ]
	orthoFile.close()
	return methAr, unmethAr

def getFPKM( fpkmFileStr ):
	'''
		returns the array or genes ordered by fpkm
	'''
	# check if exists
	checkFPKM( fpkmFileStr )
	
	fpkmFile = open( fpkmFileStr + ".txt", 'r' )
	fpkmAr = []
	fpkmValAr = []
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		fpkmAr += [lineAr[0]]
		fpkmValAr += [ lineAr[1] ]
	fpkmFile.close()
	return fpkmAr, fpkmValAr

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

def trimGenic( geneAr, fpkmAr, fpkmValAr ):
	'''
		selects numFeatures input genes from geneAr
		genes are then ordered by FPKM ranking
	'''
	
	tmpAr = []
	
	for i in range(len(geneAr)):
		try:
			fpkm = fpkmAr.index( geneAr[i][0] )
			tmpAr += [ (fpkm, geneAr[i][0], geneAr[i][1], fpkmValAr[ fpkm ] ) ]
		except ValueError:
			continue
	
	#tmpAr.sort(reverse=True)
	tmpAr.sort()
			
	# separate
	fpkmAr = []
	outAr1= []
	outAr2 = []
	for i in tmpAr:
		fpkmAr += [ i[3] ]
		outAr1 += [ i[1] ]
		outAr2 += [ i[2] ]
	return [ outAr1, outAr2 ], fpkmAr

def processBedRegions( bedFileStr, interUnmethAr, interMethAr, geneUnmethAr, geneMethAr, gffDict, numBins, numBinsStream, streamSize ):
	# out: interUnmethRegion, interMethRegion, geneUnmethName, geneMethName, rowNum, ...lots of columns...
	
	# get bedDict
	print( 'Reading {:s}...'.format( bedFileStr ) )
	bedDict, readCounts = readBed( bedFileStr )
	
	nameMatrix = []
	outMatrix = []
	numList = []
	metaAr = []
	
	metaTmpAr = [0]*numBins
	for i in range( len( interUnmethAr ) ):
		# handle interUnmeth
		# region: (chrm, start, end, numCGs, CHG methylation, CHH methylation)
		iuRegion = interUnmethAr[i]
		iuName = '{:s}:{:d}-{:d}'.format( iuRegion[0], iuRegion[1], iuRegion[2] )
		nameTmpAr = [ iuName, 'intergenic_unmeth' ]
		iuBedAr = varByRegion( bedDict, iuRegion[0], iuRegion[1], iuRegion[2], numBins )
		numTmpAr = [ x / readCounts * 1000000 for x in iuBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaTmpAr = addToArray( metaTmpAr, numTmpAr )
	metaTmpAr = [ x / len(interUnmethAr) for x in metaTmpAr ]
	metaAr += metaTmpAr
	# end for interUnmeth
	
	metaTmpAr = [0]*numBins
	for i in range( len( interMethAr ) ):
		# handle interMeth
		imRegion = interMethAr[i]
		imName = '{:s}:{:d}-{:d}'.format( imRegion[0], imRegion[1], imRegion[2] )
		imBedAr = varByRegion( bedDict, imRegion[0], imRegion[1], imRegion[2], numBins )
		nameTmpAr = [ imName, 'intergenic_meth' ]
		numTmpAr = [ x / readCounts * 1000000 for x in imBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaTmpAr = addToArray( metaTmpAr, numTmpAr )
	metaTmpAr = [ x / len(interMethAr) for x in metaTmpAr ]
	metaAr += metaTmpAr
	# end for interMeth
	
	metaTmpAr = [0]*(numBins+2*numBinsStream)
	for i in range( len( geneUnmethAr ) ):
		# handle geneUnmeth
		guName = geneUnmethAr[i]
		guBedAr = handleBedGene( bedDict, guName, gffDict, numBins, numBinsStream, streamSize )
		nameTmpAr = [ guName, 'nonGBM' ]
		numTmpAr = [ x / readCounts * 1000000 for x in guBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaTmpAr = addToArray( metaTmpAr, numTmpAr )
	metaTmpAr = [ x / len(geneUnmethAr) for x in metaTmpAr ]
	metaAr += metaTmpAr
	# end for geneUnmeth
	
	metaTmpAr = [0]*(numBins+2*numBinsStream)
	for i in range( len( geneMethAr ) ):
		# handle geneMeth
		gmName = geneMethAr[i]
		gmBedAr = handleBedGene( bedDict, gmName, gffDict, numBins, numBinsStream, streamSize )
		nameTmpAr = [ gmName, 'GBM' ]
		numTmpAr = [ x / readCounts * 1000000 for x in gmBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaTmpAr = addToArray( metaTmpAr, numTmpAr )
	metaTmpAr = [ x / len(geneMethAr) for x in metaTmpAr ]
	metaAr += metaTmpAr
	#  end for geneUnmeth
		
	return nameMatrix, outMatrix, numList, metaAr
	
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
	readCounts = 0
	
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
			readCounts += 1
	
	bedFile.close()
	return bedDict, readCounts

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

def handleBedGene( bedDict, geneName, gffDict, numBins, numBinsStream, streamSize ):
	
	geneRegion = gffDict.get( geneName )
	if geneRegion == None:
		print( 'ERROR: gene {:s} does not exist'.format( geneName ) )
	chrm = geneRegion[0]
	start = geneRegion[1]
	end = geneRegion[2]
	strand = geneRegion[3]
	
	upstreamBedAr = []
	downstreamBedAr = []
	if streamSize != 0 and numBinsStream != 0:
		upstream = start - streamSize
		downstream = end + streamSize
		upstreamBedAr = varByRegion( bedDict, chrm, upstream, start, numBinsStream )
		downstreamBedAr = varByRegion( bedDict, chrm, end, downstream, numBinsStream )
	geneBedAr = varByRegion( bedDict, chrm, start, end, numBins )
	
	if strand == '+':
		outAr = upstreamBedAr + geneBedAr + downstreamBedAr
	else:
		# upstream <- reverse downstream
		downstreamBedAr.reverse()
		# gene <- reverse gene
		geneBedAr.reverse()
		# downstream <- reverse upstream
		upstreamBedAr.reverse()
		outAr = downstreamBedAr + geneBedAr + upstreamBedAr
	return outAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	
	return oldAr

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

def writeOutput(outFileStr, nameMatrix, outMatrix, thresh, numBins, numBinsStream, info, fpkmAr):
	outFile = open( outFileStr, 'w' )
	# header: interUnmethRegion, interMethRegion, geneUnmethName, geneMethName, rowNum, 'interUnmethX' * numBins, 'interMethX' * numBins, 'geneUnmethX' * (numBinsStream*2+numBins), 'geneMethX' * (numBinsStream*2+numBins)
	# header: rowNum, gene, isGBM, FPKM, binNum
	headerAr = [ 'name', 'region', 'FPKM', 'binNum', 'value' ]
	outFile.write( info + '\n' )
	outFile.write( '\t'.join( headerAr ) + '\n' )
	
	# loop through matrices
	for i in range( len( outMatrix ) ):
		names = nameMatrix[i]
		#print( names )
		vals = outMatrix[i]
		for j in range(len(vals)):
			if vals[j] > thresh:
				vals[j] = thresh
			outStr = '{:s}\t{:s}\t{:d}\t{:.5f}\n'.format( '\t'.join( names ), fpkmAr[i], j, vals[j] )
			outFile.write( outStr )
	outFile.close()


def writeMetaOutput( outFileStr, metaMatrix, speciesIDAr, numBins, numBinsStream, info ):
	outFile = open( outFileStr, 'w' )
	# header: bin, type, species, value
	header = 'species\tregion\tbin\tvalue\n'
	outFile.write( info + '\n' + header )
	regionAr =['interUnmeth']*(numBins) + ['interMeth']*(numBins) + ['geneUnmeth']*(numBins+2*numBinsStream) + ['geneMeth']*(numBins+2*numBinsStream)
	binAr = list(range(0,numBins)) + list(range(0,numBins)) + list(range(0,numBins+2*numBinsStream)) + list(range(0,numBins+2*numBinsStream))
	# loop through species
	for i in range(len(speciesIDAr)):
		# loop through bins
		for j in range(len(regionAr)):
			outStr = '{:s}\t{:s}\t{:d}\t{:.4f}\n'.format( speciesIDAr[i], regionAr[j], binAr[j], metaMatrix[i][j] )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	
	# fpkmFileStr, orthologFileStr, numFeatures, outID, speciesIDAr, numProc
	# gffFileAr, intergenicBedAr, allcFileAr, chipBedAr, fastaFileAr, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile
	
	fpkmFileStr = None
	gffFileAr = [None, None]
	intergenicBedAr = [None, None]
	allcFileAr = [None, None]
	chipBedAr = [None, None]
	fastaFileAr = [None, None]
	orthologFileStr = None
	speciesIDAr = ['1','2']
	outID = 'out'
	numBins = NUMBINS
	numBinsStream = NUMBINSSTREAM
	regionWidth = REGIONWIDTH
	numFeatures = NUMSAMPLES
	numCGs = NUMCGS
	streamSize = STREAMSIZE
	percentile = PERCENTILE
	numProc = NUMPROC
	
	for i in range(len(argv)):
		# outID
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
		# numProc
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		# regionWidth
		elif argv[i].startswith( '-r=' ):
			try:
				regionWidth = int( argv[i][3:] )
			except ValueError:
				print( 'ERROR: intergenic region width must be integer' )
				exit()
		# numFeatures
		elif argv[i].startswith( '-n=' ):
			try:
				numFeatures = int( argv[i][3:] )
			except ValueError:
				print( 'ERROR: number of features must be integer' )
				exit()
		# numBins, numBinsStream
		elif argv[i].startswith( '-b=' ):
			try:
				numBinsAr = argv[i][3:].split(',')
				numBins = int( numBinsAr[0] )
				if len( numBinsAr ) == 2:
					numBinsStream = int( numBinsAr[1] )
				else:
					numBinsStream = int( numBinsAr[0] )
			except ValueError:
				print( 'ERROR: number of bins must be integer' )
				exit()
		# streamSize
		elif argv[i].startswith( '-s=' ):
			try:
				streamSize = int( argv[i][3:] )
			except ValueError:
				print( 'ERROR: genic upstream/downstream size must be integer' )
				exit()
		# percentile
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		# speciesIDAr[0]
		elif argv[i].startswith('-id1='):
			speciesIDAr[0] = argv[i][5:]
		# speciesIDAr[1]
		elif argv[i].startswith('-id2='):
			speciesIDAr[1] = argv[i][5:]
		# intergenicBedAr[0]
		elif argv[i].startswith('-i1='):
			intergenicBedAr[0] = argv[i][4:]
		# intergenicBedAr[1]
		elif argv[i].startswith('-i2='):
			intergenicBedAr[1] = argv[i][4:]
		# allcFileAr[0]
		elif argv[i].startswith('-a1='):
			allcFileAr[0] = argv[i][4:]
		# allcFileAr[1]
		elif argv[i].startswith('-a2='):
			allcFileAr[1] = argv[i][4:]
		# chipBedAr[0]
		elif argv[i].startswith('-b1='):
			chipBedAr[0] = argv[i][4:]
		# chipBedAr[1]
		elif argv[i].startswith('-b2='):
			chipBedAr[1] = argv[i][4:]
		# fastaFileAr[0]
		elif argv[i].startswith('-f1='):
			fastaFileAr[0] = argv[i][4:]
		# fastaFileAr[1]
		elif argv[i].startswith('-f2='):
			fastaFileAr[1] = argv[i][4:]
		# gffFileAr[0]
		elif argv[i].startswith('-g1='):
			gffFileAr[0] = argv[i][4:]
		# gffFileAr[1]
		elif argv[i].startswith('-g2='):
			gffFileAr[1] = argv[i][4:]
		# orthoFileStr
		elif argv[i].startswith( '-ortho=' ):
			orthoFileStr = argv[i][7:]
		elif argv[i].startswith( '-fpkm=' ):
			fpkmFileStr = argv[i][6:]
	# end for
	checkInputs( allcFileAr, chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr )
	processInputs( allcFileAr, chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr, intergenicBedAr, speciesIDAr, outID, numProc, regionWidth, numFeatures, numCGs, numBins, numBinsStream, streamSize, percentile )

def checkInputs( allcFileAr, chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr ):
	if orthoFileStr == None:
		print( 'ERROR: must specify otholog file. Use -ortho=<file>' )
		exit()
	elif fpkmFileStr == None:
		print( 'ERROR: must specify fpkm file for species 1. Use -fpkm=<file>' )
		exit()
	elif None in allcFileAr:
		ind = allcFileAr.index( None )
		print( 'ERROR: missing allC file for species #{:d}. Use -a{:d}=<file>'.format( ind, ind+1 ) )
		exit()
	elif None in chipBedAr:
		ind = chipBedAr.index( None )
		print( 'ERROR: missing ChIP BED file for species #{:d}. Use -b{:d}=<file>'.format( ind, ind+1 ) )
		exit()
	elif None in fastaFileAr:
		ind = fastaFileAr.index( None )
		print( 'ERROR: missing fasta file for species #{:d}. Use -f{:d}=<file>'.format( ind, ind+1 ) )
		exit()
	elif None in gffFileAr:
		ind = gffFileAr.index( None )
		print( 'ERROR: missing GFF file for species #{:d}. Use -g{:d}=<file>'.format( ind, ind+1 ) )
		exit()
	return True
		

def printHelp():
	# Usage: python3.4 chip_ortholog_compare.py [-o=outID] [-p=num_proc] [-r=intergenic_region_width] [-n=num_features] [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id] [-i1=intergenic_file_species1] [-i2=<intergenic_species2] -a1=<allc_species1> -a2=<allc_species2> -b1<chip_bed_species1> -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>
	print( 'Usage: python3.4 chip_ortholog_compare.py [-o=outID] [-p=num_proc] [-r=intergenic_region_width] [-n=num_features] [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id] [-i1=intergenic_file_species1] [-i2=<intergenic_species2] -a1=<allc_species1> -a2=<allc_species2> -b1<chip_bed_species1> -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>' )
	print( "-o=outID\toutfile identifier (default: 'out')" )
	print( '-p=num_proc\tnumber of processors to use (default: 1)' ) 
	print( '-r=intergenic_region_width\tbp width of intergenic region sections (default: 1000)' )
	print( '-n=num_features\tnumber of features to include for the heatmaps (default: 1000)' )
	print( '-b=num_bins[,num_bins_stream]\tnumber of bins for features and optionally different number of bins for upstream/downstream of genes (default 20, 0)' )
	print( '-s=gene_stream_size\tupstream/downstream size for genes (bp) (default: 0)' )
	print( '-t=percentile\tpercentile to use for heatmap correction default: (0.95)' )
	print( "-id#=species_id\tidentifier for species (default: #)" )
	print( '-i#=intergenic_file_species\tBED format file of intergenic regions for species # (default: None)' )
	print( '-a#=allc_file\tallC file for species #; file has all chromosomes' )
	print( '-b#=chip_bed_file\tBED file for species #; this is this ChIP of output heatmap' )
	print( '-f#=fasta_file\tFASTA file for species #' )
	print( '-g#=gff_file\tGFF file for species #' )
	print( '-ortho=ortholog_file\tortholog file with both species gene names and if that gene is gene body methylated' )
	print( '-fpkm=fpkm_file\tFPKM file of species 1 used to order genes for output' )

if __name__ == "__main__":
	if len(sys.argv) < 11 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
