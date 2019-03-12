import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: metaplot_expression_bed_pe.py [-p=num_proc] [-o=out_id] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-q=num_exp_groups] <gff_file> <fpkm_file> <bed_file> [bed_file]*
# create metaplot of input BED files based on expression of genes; groups genes based on number of expression groups 

NUMPROC=1
NUMBINS=20
STREAMSIZE=100
NUMEXPGROUPS=1

def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, outID, numGroups ):
	
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Reading FPKM file...' )
	fpkmAr, fpkmCount =  getFPKM( fpkmFileStr )
	print( 'LEN:', fpkmCount )
	fpkmGroupAr, fpkmCalAr = genesToGroups( fpkmAr, fpkmCount, numGroups )
	print( 'Reading GFF file...' )
	geneDict, chrmFormat = readGFF( gffFileStr )
	countGroupSt = [ '{:d}-({:d},{:.4f})'.format(i+1, len(fpkmGroupAr[i]), fpkmCalAr[i] ) for i in range(len(fpkmGroupAr)) ]
	# info
	info = "#from_script:metaplot_expression_bed_pe.py; num_bins:{:d}; num_bins_stream:{:d}; stream_size:{:d}; expression_groups(#genes,ave_fpkm):{:s}".format( numBins, numBinsStream, streamSize, ','.join( countGroupSt) )
	
	outFileStr = "{:s}_meta_expression.tsv".format( ('out' if outID == None else outID) )
	# process file -> then groups
	print( 'Begin processing {:d} files with {:d} processors...'.format( len(bedFileStrAr), numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(fpkmGroupAr, geneDict, f, numBins, numBinsStream, streamSize, chrmFormat) ) for f in bedFileStrAr ]
	# array of dictionaries -> dictionary for each sample with keys indicating expression level
	dictMatrix = [ p.get() for p in results ]
	if len(dictMatrix) != len(sampleNamesAr ):
		print( 'ERROR: not enough dictionaries for samples' )
	
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, dictMatrix, sampleNamesAr, info )
	print( 'Done.' )

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr

def getFPKM( fpkmFileStr ):
	'''
		returns the array or genes ordered by fpkm (highest first)
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
		if line.startswith( 'Potri')==False and i != -1:
			lineAr[0] = lineAr[0][:i]
		fpkmAr += [ (lineAr[0], float(lineAr[1])) ]
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
	'''
		splits genes into evently sized groups
		groups are ordered with highest FPKM values first
	'''
	outAr = []
	calcAr = []
	numPerGroup = int( math.ceil( fpkmCount / numGroups ) )
	tmpAr = []
	tmpCal = 0
	# loop through genes by fpkm value
	#print( fpkmAr )
	for i in range(len(fpkmAr)):
		# enough for group
		if len(tmpAr) >= numPerGroup:
			#print( 'tmp', tmpAr )
			outAr.append( tmpAr )
			calcAr += [ tmpCal / len(tmpAr) ]
			tmpAr = [ fpkmAr[i][0] ]
			tmpCal = fpkmAr[i][1]
			print( i, fpkmAr[i][0], fpkmAr[i][1] )
		else:
			tmpAr += [ fpkmAr[i][0] ]
			tmpCal += fpkmAr[i][1]
	# add last group
	outAr.append( tmpAr )
	calcAr += [ tmpCal / len(tmpAr) ]
	return outAr, calcAr

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
			#print( name )
			# put into geneDict
			gffDict[name] = (chrm, start, end, strand)
	
	gffFile.close()
	return gffDict, chrmFormat

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return notesStr[adIndex:endIndex+adIndex]

def processFile( fpkmGroupAr, gffDict, fileStr, numBins, numBinsStream, streamSize, chrmFormat ):
	print( 'Reading BED file {:s}...'.format( fileStr ) )
	bedDict, countReads = readBed( fileStr )
	#print( 'Found {:d} reads.'.format( countReads ) )
	print( 'Processing {:s}...'.format( fileStr ) )
	outDict = {}
	n = len(fpkmGroupAr)
	for i in range(n):
		#print( 'Group {:d} - {:d}'.format( i+1, len(fpkmGroupAr[i] )))
		outDict [i+1] = processBed( fpkmGroupAr[i], gffDict, bedDict, countReads, numBins, numBinsStream, streamSize )
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
		#print( gene )
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
	outAr = [ float(x)/(countReads*geneCount)*1000000000 for x in outAr ]
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
	varAr = adjustArrayLength( binWidth, varAr )
	return varAr

def adjustArrayLength( length, inAr ):
	outAr = [0] * len( inAr )
	for i in range( len(inAr) ):
		outAr[i] = float( inAr[i] ) / float( length )
	return outAr

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
	outFile.write( '{:s}\nsample\texpression\tbin\tvalue\n'.format( info ) )
	
	# for each bed file
	for n in range(len(sampleNamesAr )):
		# for each expression level
		for ex in sorted( dictMatrix[n].keys() ):
			tmpAr = dictMatrix[n][ex]	# array of bins w/ values
			outAr = [ '{:s}\t{:d}\t{:d}\t{:f}'.format( sampleNamesAr[n], ex, i+1, tmpAr[i] ) for i in range(len(tmpAr) ) ]
			outFile.write( '\n'.join( outAr ) + '\n' )
		# end for ex
	# end for bed
	outFile.close()
	
		
def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	numGroups = NUMEXPGROUPS
	outID = None
	startInd = 0

	for i in range(min(5,len(argv)-3)):
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
		elif argv[i].startswith( '-q=' ):
			try:
				numGroups = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of expression groups must be integer' )
				exit()
	# end for
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	
	bedFileStrAr = []
	for j in range(startInd+2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, outID, numGroups )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 metaplot_expression_bed_pe.py [-p=num_proc] [-o=out_id] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-q=num_exp_groups] <gff_file> <fpkm_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
