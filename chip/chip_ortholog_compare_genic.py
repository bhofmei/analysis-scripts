import sys, math, glob, multiprocessing, subprocess, os, random

# Usage: python3.4 chip_ortholog_compare.py [-o=outID] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id] -b1<chip_bed_species1> -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>
# Using methylation information and gene orthologs creates heatmaps showing ChIP enrichment for two samples; includes intergenic and genic regions

# allC_files chip_bed_files ortholog_file gff_files intergenic_bed_files
# fasta_files


NUMBINS=20
NUMBINSSTREAM=0
STREAMSIZE=0
PERCENTILE=0.95
NUMPROC=1

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

def processInputs( chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr, speciesIDAr, outID, numProc, numBins, numBinsStream, streamSize, percentile ):
	# fpkmFileStr, orthologFileStr, numFeatures, outID, speciesIDAr, numProc
	# gffFileAr, intergenicBedAr, allcFileAr, chipBedAr, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile
	
	gffStr = [ os.path.basename( x ) for x in gffFileAr ]
	chipBedStr = [ os.path.basename( x ) for x in chipBedAr ]
	# write meta information
	print( numBins, numBinsStream, streamSize, os.path.basename(orthoFileStr), gffStr, chipBedStr, os.path.basename(fpkmFileStr) )
	info = '#from_script:chip_ortholog_compare_genic.py numBins:{:d}; numBinsStream:{:d}; streamSize:{:d}; ortho_file:{:s}; gff_files:{:s}; chip_bed_files:{:s}; fpkm_file:{:s}'.format( numBins, numBinsStream, streamSize, os.path.basename(orthoFileStr), ','.join(gffStr), ','.join(chipBedStr), os.path.basename(fpkmFileStr) )
	
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
	# gffFileStr, chipBedStr, outFileStr, unmethAr, methAr, speciesID, numBins, numBinsStream, streamSize, percentile
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing species with {:d} processes'.format( numProc ) )
	results = [ pool.apply_async( processSpecies, args=( gffFileAr[i], chipBedAr[i], fastaFileAr[i], outFileAr[i], unmethArs[i], methArs[i], speciesIDAr[i], numBins, numBinsStream, streamSize, percentile, fpkmMethAr, fpkmUnmethAr, info ) ) for i in range(len(speciesIDAr)) ]
	metaMatrix = [ p.get() for p in results ]
	outMetaStr = '{:s}_meta.tsv'.format( outID )
	print( 'Writing metaplot output to {:s}...'.format( outMetaStr ) )
	writeMetaOutput( outMetaStr, metaMatrix, speciesIDAr, numBins, numBinsStream, info )
	print( 'Done.' )

def processSpecies( gffFileStr, chipBedStr, fastaFileStr, outFileStr, unmethAr, methAr, speciesID, numBins, numBinsStream, streamSize, percentile, fpkmMethAr, fpkmUnmethAr, info ):

	# gff dict
	print( 'Reading {:s}...'.format( gffFileStr ) )
	gffDict = readGFF( gffFileStr )
	
	# process regions
	print( 'Processing all regions for {:s}...'.format( speciesID ) )
	nameMatrix, outMatrix, numList, metaAr = processBedRegions( chipBedStr, unmethAr, methAr, gffDict, numBins, numBinsStream, streamSize )
	# threshold
	thresh = calculatePercentile( percentile, numList )
	if numBinsStream != 0 and streamSize == 0:
		numBinsStream = 0
	fpkmAr = fpkmUnmethAr + fpkmMethAr
	# write output
	print( 'Writing species {:s} output to {:s}'.format( speciesID, outFileStr ) )
	#info = '#numFeatures:{:d}_regionWidth:{:d}_numCGs:{:d}_numBins:{:d}_numBinsStream:{:d}_streamSize:{:d}_percentile:{:.3f}'.format( numFeatures, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile )
	
	writeOutput( outFileStr, nameMatrix, outMatrix, thresh, numBins, numBinsStream, info, fpkmAr)
	return metaAr

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
	isFirst = False
	
	orthoFile = open( orthologFileStr, 'r' )
	for line in orthoFile:
		if isFirst:
			isFirst = False
			continue
		lineAr = line.rstrip().split(  )
		# (0) Athaliana_gene (1) gBM/UM (2) Esalsugineum_gene (3) gBM/UM
		atGene = lineAr[0][:-2]
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

def processBedRegions( bedFileStr, geneUnmethAr, geneMethAr, gffDict, numBins, numBinsStream, streamSize ):
	# out: geneUnmethName, geneMethName, rowNum, ...lots of columns...
	
	# get bedDict
	print( 'Reading {:s}...'.format( bedFileStr ) )
	bedDict, readCounts = readBed( bedFileStr )
	
	nameMatrix = []
	outMatrix = []
	numList = []
	metaArU = [0] * (numBins + 2*numBinsStream )
	metaArG = [0] * (numBins + 2*numBinsStream )
	
	
	# loop through unmethylated genes
	for i in range( len( geneUnmethAr ) ):
		# handle geneUnmeth
		guName = geneUnmethAr[i]
		guBedAr = handleBedGene( bedDict, guName, gffDict, numBins, numBinsStream, streamSize )
		
		nameTmpAr = [ guName, 'nonGBM' ]
		# adjust by readCounts
		numTmpAr = [ x / readCounts * 1000000 for x in guBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaArU = addToArray( metaArU, numTmpAr )
	#end for unmethylated
	
	
	for i in range( len( geneMethAr ) ):
		# handle geneMeth
		gmName = geneMethAr[i]
		gmBedAr = handleBedGene( bedDict, gmName, gffDict, numBins, numBinsStream, streamSize )
		nameTmpAr = [ gmName, 'GBM' ]
		# adjust by readCounts
		numTmpAr = [ x / readCounts * 1000000 for x in gmBedAr ]
		nameMatrix.append( nameTmpAr )
		outMatrix.append( numTmpAr )
		numList += numTmpAr
		metaArG = addToArray( metaArG, numTmpAr )
	#end for unmethylated
	metaAr = metaArU + metaArG
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
	headerAr = [ 'geneName', 'gBM', 'FPKM', 'binNum', 'value' ]
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
	regionAr = ['geneUnmeth']*(numBins+2*numBinsStream) + ['geneMeth']*(numBins+2*numBinsStream)
	# loop through species
	for i in range(len(speciesIDAr)):
		# loop through bins
		for j in range(len(regionAr)):
			outStr = '{:s}\t{:s}\t{:d}\t{:.4f}\n'.format( speciesIDAr[i], regionAr[j], j, metaMatrix[i][j] )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	
	# fpkmFileStr, orthologFileStr, numFeatures, outID, speciesIDAr, numProc
	# gffFileAr, intergenicBedAr, allcFileAr, chipBedAr, fastaFileAr, regionWidth, numCGs, numBins, numBinsStream, streamSize, percentile
	
	fpkmFileStr = None
	gffFileAr = [None, None]
	chipBedAr = [None, None]
	fastaFileAr = [None, None]
	orthologFileStr = None
	speciesIDAr = ['1','2']
	outID = 'out'
	numBins = NUMBINS
	numBinsStream = NUMBINSSTREAM
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
	checkInputs( chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr )
	processInputs( chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr, speciesIDAr, outID, numProc, numBins, numBinsStream, streamSize, percentile )

def checkInputs( chipBedAr, fastaFileAr, gffFileAr, orthoFileStr, fpkmFileStr ):
	if orthoFileStr == None:
		print( 'ERROR: must specify otholog file. Use -ortho=<file>' )
		exit()
	elif fpkmFileStr == None:
		print( 'ERROR: must specify fpkm file for species 1. Use -fpkm=<file>' )
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
	# Usage: python3.4 chip_ortholog_compare.py [-o=outID] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id] -b1<chip_bed_species1> -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>
	print( 'Usage: python3.4 chip_ortholog_compare.py [-o=outID] [-p=num_proc]  [-b=num_bins[,num_bins_stream]] [-s=gene_stream_size] [-t=percentile] [-id1=species1_id] [-id2=species2_id]  -b2=<chip_bed_species2> -f1=<fasta_species1> -f2=<fasta_species2> -g1=<gff_species1> -g2=<gff_species2>  -ortho=<ortholog_file> -fpkm=<fpkm_species1>' )
	print( "-o=outID\toutfile identifier (default: 'out')" )
	print( '-p=num_proc\tnumber of processors to use (default: 1)' ) 
	print( '-b=num_bins[,num_bins_stream]\tnumber of bins for features and optionally different number of bins for upstream/downstream of genes (default 20, 0)' )
	print( '-s=gene_stream_size\tupstream/downstream size for genes (bp) (default: 0)' )
	print( '-t=percentile\tpercentile to use for heatmap correction default: (0.95)' )
	print( "-id#=species_id\tidentifier for species (default: #)" )
	print( '-b#=chip_bed_file\tBED file for species #; this is this ChIP of output heatmap' )
	print( '-f#=fasta_file\tFASTA file for species #' )
	print( '-g#=gff_file\tGFF file for species #' )
	print( '-ortho=ortholog_file\tortholog file with both species gene names and if that gene is gene body methylated' )
	print( '-fpkm=fpkm_file\tFPKM file of species 1 used to order genes for output' )

if __name__ == "__main__":
	if len(sys.argv) < 9 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
