import sys, math, glob, multiprocessing, subprocess, os, bisect

# Usage: python3.4 metaplot_ortho_expression_bed_pe.py [-1|-2] [-a|-r|-g|-n] [-p=num_proc] [-o=out_id] [-f=ortholog_file] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-q=num_exp_groups] <gff_file> <fpkm_file> <bed_file> [bed_file]*

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000

ORTHO_FILE='/Users/bhofmeister/Documents/Research/collaborations/cmt3/arabidopsis_compare/Esalsugineum_Athaliana_bin.out'
def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, outID, numProc, numBins, numBinsStream, streamSize, genesIncluded, orthoFile, numGroups, speciesInt ):
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nGenes included: {:s}\nNumber of expression groups: {:d}\nSpecies: {:d}\nSamples included: {:s}\n'.format( numBins, numBinsStream, streamSize, genesIncluded, numGroups, speciesInt, ', '.join(sampleNamesAr) ) )
	if genesIncluded == 'all':
		subsetAr = None
	else:
		print( 'Reading A. thaliana - E. salsugineum ortholog file...' )
		subsetAr = readGBMPairs( genesIncluded, orthoFile, speciesInt )
	print( 'Reading FPKM file...' )
	fpkmAr, fpkmCount =  getFPKM( fpkmFileStr, subsetAr )
	print( 'LEN:', fpkmCount )
	fpkmGroupAr = genesToGroups( fpkmAr, fpkmCount, numGroups )
	infoStr = '#'
	for i in range(len(fpkmGroupAr)):
		infoStr += '{:d}-{:d}..'.format(len(fpkmGroupAr)-i,len(fpkmGroupAr[i]))
	infoStr += 'total-{:d}'.format( int(fpkmCount))
	#print (fpkmGroupAr)
	print( 'Reading GFF file...' )
	geneDict, chrmFormat = readGFF( gffFileStr )
	
	if outID == '' and genesIncluded == 'all':
		outFileStr = 'meta_{:d}_{:d}.csv'.format(numBins,streamSize)
	elif outID == '':
		outFileStr = 'meta_{:d}_{:d}_{:s}.csv'.format(numBins,streamSize,genesIncluded)
	elif genesIncluded == 'all':
		outFileStr = 'meta_{:d}_{:d}_{:s}.csv'.format(numBins,streamSize,outID)
	else:
		outFileStr = 'meta_{:d}_{:d}_{:s}_{:s}.csv'.format(numBins,streamSize,outID,genesIncluded)
	
	print( 'Begin processing files with {:d} processors...'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(fpkmGroupAr, geneDict, f, numBins, numBinsStream, streamSize, chrmFormat) ) for f in bedFileStrAr ]
	# array of dictionaries -> dictionary for each sample with keys indicating expression level
	dictMatrix = [ p.get() for p in results ]
	if len(dictMatrix) != len(sampleNamesAr ):
		print( 'ERROR: not enough dictionaries for samples' )
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, dictMatrix, sampleNamesAr, infoStr )
	print( 'Done.' )

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr
	
def readGBMPairs( genesIncluded, orthoFile, speciesInt ):
	
	outAr = []
	inFile = open( orthoFile, 'r' )
	ind = speciesInt * 2 + -2
	for line in inFile:
		lineAr = line.rstrip().split()
		# remove .# from athaliana cds
		i = lineAr[0].rfind('.')
		if i != -1:
			lineAr[0] = lineAr[0][:i]
		# (0) athaliana cds (1) UM/CG (2) eutrema cds (3) UM/CG
		# ortholog
		if genesIncluded == 'orthologs':
			bisect.insort( outAr, lineAr[ind] )
		# gbm
		elif genesIncluded == 'gBM_orthologs' and lineAr[1] == 'CG':
			bisect.insort( outAr, lineAr[ind] )
		# non-gbm
		elif genesIncluded == 'non_gBM_orthologs' and lineAr[1] == 'UM':
			bisect.insort( outAr, lineAr[ind] )
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
	fpkmCount = 0
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		i = lineAr[0].rfind('.')
		if i != -1:
			lineAr[0] = lineAr[0][:i]

		# need gene
		if subsetAr == None or bisectIndex( subsetAr, lineAr[0] ) != None:
			fpkmAr += [ lineAr[0] ]
			fpkmCount += 1.0
	fpkmFile.close()
	return fpkmAr, fpkmCount

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

def genesToGroups( fpkmAr, fpkmCount, numGroups ):
	
	outAr = []
	numPerGroup = int( math.ceil( fpkmCount / numGroups ) )
	tmpAr = []
	# loop through genes by fpkm value
	#print( fpkmAr )
	for i in range(len(fpkmAr)):
		# enough for group
		if len(tmpAr) >= numPerGroup:
			#print( 'tmp', tmpAr )
			outAr.append( tmpAr )
			tmpAr = [ fpkmAr[i] ]
		else:
			tmpAr += [ fpkmAr[i] ]
	# add last group
	outAr.append( tmpAr )
	return outAr

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
	
def processFile( fpkmGroupAr, gffDict, fileStr, numBins, numBinsStream, streamSize, chrmFormat ):
	print( 'Reading BED file {:s}...'.format( fileStr ) )
	bedDict, countReads = readBed( fileStr )
	print( 'Processing {:s}...'.format( fileStr ) )
	outDict = {}
	n = len(fpkmGroupAr)
	for i in range(n):
		print( 'Group {:d} - {:d}'.format( n-i, len(fpkmGroupAr[i] )))
		outDict [n-i] = processBed( fpkmGroupAr[i], gffDict, bedDict, countReads, numBins, numBinsStream, streamSize )
	print( 'Finished processing', fileStr )
	return outDict

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
	countReads = 0
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		lineAr =line.rstrip().split('\t')
		try:
			curScaf = lineAr[0]
			pos = int( lineAr[1] ) + 1
		except ValueError:
			pass
		else:
			countReads += 1
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
	return bedDict, countReads

	
def processBed( geneAr, gffDict, bedDict, countReads, numBins, numBinsStream, streamSize ):
	# repeatAr organized [(scaffold, start, end, strand), ...]
	
	totalGeneAr = [0] * numBins
	totalUpAr = []
	totalDownAr = []
	geneCount = 0
	
	if streamSize != 0:
		totalUpAr = [0] * numBinsStream
		totalDownAr = [0] * numBinsStream
	
	# loop by gene length
	for gene in geneAr:
		info = gffDict.get( gene )
		# not actually gene - don't count it
		#print( region )
		if info == None:
			continue
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		geneCount += 1
		outAr = []
		if streamSize != 0:
			upstream = start - streamSize
			downstream = end + streamSize
			# upstream
			curUpVarAr = varByRegion(bedDict, chrm, upstream, start, numBinsStream)
			# downstream
			curDownVarAr = varByRegion(bedDict, chrm, end, downstream, numBinsStream)
		
		# gene region
		curGeneVarAr = varByRegion(bedDict, chrm, start, end, numBins)
		#print( curGeneVarAr )
		# forward strand - all arrays are okay
		if strand == "+":
			if streamSize != 0:
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		else:
			if streamSize != 0:
				# upstream <- reverse downstream arrays
				curDownVarAr.reverse()
				# gene <- reverse gene arrays
				curGeneVarAr.reverse()
				# downstream <- reverse upstream arrays
				curUpVarAr.reverse()
				#outAr = curDownVarAr + curRepeatVarAr + curUpVarAr
				totalUpAr = addToArray( totalUpAr, curUpVarAr )
				totalDownAr = addToArray( totalDownAr, curDownVarAr )
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
			else:
				curGeneVarAr.reverse()
				totalGeneAr = addToArray( totalGeneAr, curGeneVarAr )
		# end if-else strand
	# end for gene
	outAr = totalUpAr + totalGeneAr + totalDownAr
	#print( outAr )
	outAr = [ float(x)/(countReads*geneCount)*1000000 for x in outAr ]
	return outAr

def varByRegion(bedDict, chrm, start, end, numBins):
	''' 
		takes in the variant dictionary generated for a gene and the start and 
		end positions of that gene. 
		Note the dictionary is set up so the key is the position
		on the scaffold and the value is the frequency of that variant
		returns one array with number of variants separated by bin
	'''
	binWidth = int(math.floor ( (end - start + 1) / numBins ))
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
				#print( 'ERROR with', chrm )
				pass
	
	# handle the last "catch-all" bin
	for key in range( ( (numBins-1)*binWidth+start ), ( end ) ):
		try:
			dictEntry = bedDict.get(chrm).get(key)
			# only want to include sites we have info for
			if dictEntry != None:
				varAr[-1] += dictEntry
		except AttributeError:
			#print( 'ERROR with', chrm )
			pass
	return varAr

def addToArray(oldAr, currAr):
	if len(oldAr) != len(currAr):
		return -1
	else:
		for i in range( len(oldAr) ):
			oldAr[i] += currAr[i]
	
	return oldAr
	
def writeOutput( outFileStr, dictMatrix, sampleNamesAr, info ):
	outFile = open( outFileStr, 'w' )
	
	# header
	outFile.write( 'sample,expression,bin,value\n{:s}\n'.format( info ) )
	
	# for each bed file
	for n in range(len(sampleNamesAr )):
		# for each expression level
		for ex in sorted( dictMatrix[n].keys(),reverse=True ):
			tmpAr = dictMatrix[n][ex]	# array of bins w/ values
			outAr = [ '{:s},{:d},{:d},{:f}'.format( sampleNamesAr[n], ex, i+1, tmpAr[i] ) for i in range(len(tmpAr) ) ]
			outFile.write( '\n'.join( outAr ) + '\n' )
		# end for ex
	# end for bed
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	genesIncluded = 'all'
	numGroups = 1
	speciesInt = 0
	outID = ''
	startInd = 0
	orthoFile = ORTHO_FILE
	
	for i in range( min(10, len(argv)-3 ) ):
		# number of processors
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
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
		elif argv[i].startswith( '-q=' ):
			try:
				numGroups = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of expression groups must be integer' )
				exit()
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
	elif speciesInt == 0:
		speciesInt = 2
	
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, outID, numProc, numBins, numBinsStream, streamSize, genesIncluded, orthoFile, numGroups, speciesInt )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 metaplot_ortho_expression_bed_pe.py [-1|-2] [-a|-r|-g|-n] [-p=num_proc] [-o=out_id] [-f=ortholog_file] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-q=num_exp_groups] <gff_file> <fpkm_file> <bed_file> [bed_file]*\n-a=all genes\n-r=all orthologous genes\n-g=all gene body methylated genes\n-n=all non-gene body methylated genes\n-1/2 species that fpkm belong to")
	else:
		parseInputs( sys.argv[1:] )
