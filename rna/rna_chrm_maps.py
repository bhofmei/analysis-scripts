import sys, math, glob, multiprocessing, subprocess, os

BINSIZE=100
THRESHOLD=2
TICKVAL=0
PERCENTILE=0.95

# Usage: python3.4 rna_chrm_maps.py [-ticks] [-o=out_id] [-c=chrms] [-b=bin_size] [-t=threshold] [-p=percentile] [-i=include_file] [-e=exclude_file] <gff_file> <fpkm_file> [wig_file]*
# note: threshold is log2 threshold

def processInputs( gffFileStr, fpkmFileStr, wigFileAr, outID, chrmList, binSize, threshold, percentile, useTicks, includeFileStr, excludeFileStr ):

	info = "#from_script:rna_chrm_maps.py;bin_size:{:d}".format( binSize )
	# read gff file -> get gene coordinates and strand
	infoGenes = info + 'gff_file:{:s}'.format( os.path.basename( gffFileStr ) )
	print('reading gff...')
	chrmDict, gffDict, geneAr = readGFF( gffFileStr )
	# geneAr contains genes that should be looked for in FPKM diagrams
	print(len(geneAr))
	if chrmList == None:
		chrmList = list(chrmDict.keys())
	else:
		infoGenes += ';chrm_list:{:s}'.format( ','.join( chrmList ) )
	print(chrmDict)
	print('processing genes')
	geneValDict= processGenes( gffDict, chrmDict, chrmList, binSize )
	# write file and delete dictionaries
	outGeneStr = outID + '_genes.tsv'
	print('writing gene information to {:s}...'.format(outGeneStr))
	writeGenes(outGeneStr, geneValDict, infoGenes)
	del gffDict, geneValDict
	
	infoFPKM = info + ';fpkm_file:{:s};threshold:{:.1f};use_ticks:{:s}'.format( os.path.basename( fpkmFileStr ) , threshold, str(useTicks))
	# clean gene list
	if includeFileStr != None:
		infoFPKM += ';include_gene_file:{:s}'.format( os.path.basename( includeFileStr ) )
		geneAr = includeGenes( geneAr, includeFileStr )
	elif excludeFileStr != None:
		infoFPKM += ';exclude_gene_file:{:s}'.format( os.path.basename( excludeFileStr) )
		geneAR = excludeGenes( geneAr, excludeFileStr )
	
	print('reading fpkm file...')
	sampleNameAr, fpkmDict = readFPKM( fpkmFileStr, geneAr )
	# handle genes and fpkm values
	print('processing FPKM' )
	fpkmValDict = processFPKM( fpkmDict, chrmDict, chrmList, binSize, sampleNameAr )
	#print( fpkmValDict )
	print( 'cleaning FPKM')
	if useTicks:
		fpkmCleanDict = cleanFPKMTicks( fpkmValDict, threshold )
	else:
		fpkmCleanDict = cleanFPKM( fpkmValDict, threshold )
	outFPKMStr = outID + '_fpkm'+('_ticks' if useTicks else '') +  '.tsv'
	print('writing fpkm information to {:s}...'.format(outFPKMStr))
	writeFPKM( outFPKMStr, fpkmCleanDict, sampleNameAr, infoFPKM )
	
	# handle bed files
	if len(wigFileAr) > 0:
		wigNames = [ getSampleName(x) for x in wigFileAr ]
		wigFiles = [os.path.basename( x ) for x in wigFileAr ]
		infoWig = info + ';wig_files:{:s};percentile:{:.1f}'.format( ','.join( wigFiles), percentile )
		print('reading wig files...')
		wigDictAr = [ readWig( x, chrmDict, chrmList, binSize ) for x in wigFileAr ]
		wigThres = [ correctWig(w, percentile) for w in wigDictAr ]
		outWigStr =  outID + '_{:d}_wig.tsv'.format(int(percentile*100))
		print('writing wig information to {:s}...'.format(outWigStr))
		writeWig( outWigStr, wigDictAr, wigNames, wigThres, infoWig )
	print( 'Done.' )
	
	
def readGFF( gffFileStr ):
	gffFile = open( gffFileStr, 'r' )
	chrmDict = {}
	geneAr = []
	gffDict = {}
	chrmSet = set()
	chrmDictGene = {}
	
	for line in gffFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		if lineAr[2] == "chromosome":
			chrmDict[ lineAr[0] ] = int(lineAr[4] )
		elif lineAr[2] == "gene":
			chrm = lineAr[0]
			start = int(lineAr[3])
			end = int(lineAr[4])
			name = getGeneName( lineAr[8] )
			if gffDict.get(chrm) == None:
				gffDict[chrm] = []
			gffDict[chrm] += [( start, end,lineAr[6] )]
			# in case chromosomes aren't part of GFF
			if chrm not in chrmSet:
				chrmSet.add(chrm)
				chrmDictGene[chrm] = end
			if end > chrmDictGene[chrm]:
				chrmDictGene[chrm] = end
		elif lineAr[2] == "mRNA":
			rName = getRNAName( lineAr[8] )
			geneAr += [ rName ]
	# end for
	if chrmDict == {}:
		chrmDict = chrmDictGene
	gffFile.close()
	return chrmDict, gffDict, geneAr

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def getRNAName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find('-')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def includeGenes( geneAr, includeFileStr ):
	tmpAr = []
	includeFile = open( includeFileStr, 'r' )
	
	for line in includeFile:
		gene = line.strip().rstrip()
		tmpAr += [ gene ]
	# for line
	includeFile.close()
	
	outAr = []
	for gene in geneAr:
		if gene in tmpAr:
			outAr += [ gene ]
	print( 'num genes total: {:d}, num_genes_included: {:d}, num_genes_left: {:d}'.format( len(geneAr), len(tmpAr), len(outAr) ) )
	return outAr

def excludeGenes( geneAr, excludeFileStr ):
	tmpAr = []
	excludeFile = open( excludeFileStr, 'r' )
	
	for line in excludeFile:
		gene = line.strip().rstrip()
		tmpAr += [ gene ]
	# for line
	excludeFile.close()
	
	outAr = []
	for gene in geneAr:
		if gene not in tmpAr:
			outAr += [ gene ]
	print( 'num genes total: {:d}, num_genes_excluded: {:d}, num_genes_left: {:d}'.format( len(geneAr), len(tmpAr), len(outAr) ) )
	return outAr

def readFPKM( fpkmFileStr, geneAr ):
	sampleNameAr = []
	indexAr = []
	locusIndex = 5
	geneIndex = 4
	isHeader = True
	fpkmDict = {}
	
	fpkmFile = open( fpkmFileStr, 'r' )
	for line in fpkmFile:
		lineAr = line.rstrip().split('\t')
		if isHeader:
			# header
			for j in range(len(lineAr)):
				if lineAr[j] == "locus":
					locusIndex = j
				elif lineAr[j].endswith('_FPKM'):
					indexAr += [j]
					sampleNameAr += [ lineAr[j].replace('_FPKM','')]
				if lineAr[j] == "gene_short_name":
					geneIndex = j
			isHeader = False
		else:
			geneName = lineAr[geneIndex]
			if geneName not in geneAr:
				continue
			locusName = lineAr[locusIndex]
			fpkmValAr = [0]*len(indexAr)
			for k in range(len(indexAr)):
				if lineAr[indexAr[k]+3]=="OK":
					fpkmValAr[k] +=  float( lineAr[indexAr[k]] ) 
			# add to dictionary
			fpkmDict[locusName] = fpkmValAr
	# end for
	fpkmFile.close()
	return sampleNameAr, fpkmDict

def getSampleName( fileStr ):
	# bed file
	leftIndex = fileStr.rfind('/')
	rightIndex = fileStr.rfind('.')
	sampleName = fileStr[leftIndex+1:rightIndex]
	return sampleName

def processGenes( gffDict, chrmDict, chrmList, binSize ):
	outDict = {}
	# loop through chromosomes
	for chrm in chrmList:
		numBins = int(math.ceil(chrmDict[chrm]/float(binSize)))
		outDict[chrm] = [0] * numBins
		# loop through genes
		for gene in gffDict[chrm]:
			start = gene[0]
			end = gene[1]
			strand = gene[2]
			bStart = start // binSize
			bEnd = end // binSize
			for bin in range(bStart, bEnd+1):
				outDict[chrm][bin] += ( 1 if strand == '+' else -1 )
	# end for chrm
	return outDict
	
def processFPKM( fpkmDict, chrmDict, chrmList, binSize, sampleNamesAr ):
	fpkmValDict = {}
	for chrm in chrmList:
		numBins = int(math.ceil(chrmDict[chrm]/float(binSize)))
		fpkmValDict[chrm] = [ [0]*numBins for i in range(len(sampleNamesAr)) ]
	# loop through genes
	for locus in fpkmDict.keys():
		#print(locus)
		l = locus.split(':')
		if len(l) != 2:
			continue
		chrm = l[0]
		if chrm not in chrmList:
			continue
		r = l[1]
		start, end = r.split('-')
		bStart = int(start) // binSize
		bEnd = int(end) // binSize
		# fpkm dictionary
		fpkmAr = fpkmDict[locus]
		#print( fpkmAr )
		# check at least one value > 1
		t = [ (1 if x < 1 else 0) for x in fpkmAr ]
		ts = sum(t)
		if ts == 4:
			y = 0
		else:
			y = 1
		# loop through sample comparison
		for j in range(len(fpkmAr)):
			val = computeChange(fpkmAr[0],fpkmAr[j])
			for bin in range(bStart, bEnd + 1):
				fpkmValDict[chrm][j][bin] = max(y * val, fpkmValDict[chrm][j][bin] )
	# end for genes
	return fpkmValDict

def computeChange( value1, value2 ):
	# use -inf and inf
	# correct for which is bigger
	if value1 == 0 and value2 == 0:
		return 0
	elif value1 == 0:
		return float('inf')
	elif value2 == 0:
		return float('-inf')
	return math.log2(value2) - math.log2(value1)
		
def readWig( fileStr, chrmDict, chrmList, binSize ):
	wigDict = {}
	wigFile = open( fileStr, 'r' )
	chrmAr = []
	chrm = None
	for line in wigFile:
		if line.startswith( 'track' ):
			pass
		# new chrm
		elif line.startswith( 'variableStep' ):
			# add previous
			if chrm != None:
				wigDict[chrm] = chrmAr
			lineAr = line.rstrip().split(' ')
			chrmT = lineAr[1].replace('chrom=','')
			if chrmDict.get(chrmT) == None or chrmT not in chrmList:
				chrm = None
			else:
				chrm = chrmT
				numBins = int(math.ceil(chrmDict[chrm]/float(binSize)))
				#print(numBins)
				chrmAr = [0] * numBins
		elif chrm != None:
			lineAr = line.rstrip().split(' ')
			bin = int( lineAr[0] ) // binSize
			#print( '--', bin )
			value = float(lineAr[1])
			chrmAr[bin] += value
	# end for
	if chrm != None:
		wigDict[chrm] = chrmAr
	wigFile.close()
	return wigDict

def correctWig( wigDict, percentile ):
	numList = []
	# loop through chromosomes and collect numbers
	for chrm in wigDict.keys():
		numList += wigDict[chrm]
	if percentile == 1:
		return max(numList)
	numList.sort()
	ind = math.ceil(percentile * len( numList ) - 1)
	try:
		return numList[ind]
	except IndexError:
		return numList[-1]

def writeGenes( outFileStr, outDict, info ):
	outFile = open( outFileStr, 'w' )
	header = info + '\n' + "#chrm\tpos\tvalue\n"
	outFile.write(header)
	# loop through chromosomes
	for chrm in sorted(outDict.keys()):
		# loop through bins
		binAr = outDict[chrm]
		for j in range(len(binAr)):
			if binAr[j] > 1:
				binAr[j] = 1
			elif binAr[j] < -1:
				binAr[j] = 1
			outFile.write( '{:s}\t{:d}\t{:d}\n'.format( chrm, j, binAr[j] ) )
	# end for chrm
	outFile.close()

def cleanFPKMTicks( outDict, threshold ):
	newDict = {}
	# loop through chms
	for chrm in sorted(outDict.keys()):
		newDict[chrm] = [ None ] *len(outDict[chrm])
		# loop through samples
		for i in range(len(outDict[chrm])):
			arr = outDict[chrm][i]
			newAr = [0] * len(arr)
			# loop through bins:
			for j in range(len(arr)):
				# check upper threshold
				if arr[j] >= threshold:
					#print( j, max(0,j-2), min(len(arr),j+3) )
					for k in range( max(0,j-TICKVAL), min(len(arr), j+TICKVAL+1)):
						newAr[k] = threshold
			# end for j
			newDict[chrm][i] = newAr
		# end for i
	# end for chrm
	return newDict

def cleanFPKM( outDict, threshold ):
	newDict = {}
	# loop through chms
	for chrm in sorted(outDict.keys()):
		newDict[chrm] = [ None ]*len(outDict[chrm])
		# loop through samples
		for i in range(len(outDict[chrm])):
			arr = outDict[chrm][i]
			# loop through bins:
			for j in range(len(arr)):
				# check upper threshold
				if arr[j] > threshold:
					arr[j] = threshold
				elif arr[j] < -1 * threshold:
					arr[j] = -1 * threshold
			# end for j
			newDict[chrm][i] = arr
		# end for i
	# end for chrm
	return newDict

def writeFPKM( outFileStr, outDict, sampleAr, info ):
	outFile = open( outFileStr, 'w' )
	header =info + '\n'+ '#chrm\tsample\tpos\tvalue\n'
	outFile.write(header)
	
	# loop through chromosomes
	for chrm in sorted(outDict.keys()):
		# loop through samples
		for i in range(len(sampleAr)):
			tmpAr = outDict[chrm][i]
			#print(tmpAr)
			# loop through bins
			for j in range(len(tmpAr)):
				outFile.write('{:s}\t{:s}\t{:d}\t{:.4f}\n'.format( chrm, sampleAr[i], j, tmpAr[j] ))
			# end for j
		# end for i
	# end for chrm
	outFile.close()

def writeWig( outFileStr, wigDictAr, wigNameAr, wigThres, info ):
	outFile = open( outFileStr, 'w' )
	header = info + '\n' + '#wig\tchrm\tpos\tvalue\n'
	outFile.write(header)
	#print(wigDictAr)
	# loop through samples
	for i in range(len(wigDictAr)):
		# loop through chrms
		for chrm in sorted(wigDictAr[i].keys()):
			tmpAr = wigDictAr[i][chrm]
			# loop through bins
			for j in range(len(tmpAr)):
				outFile.write('{:s}\t{:s}\t{:d}\t{:.4f}\n'.format( wigNameAr[i], chrm, j, ( wigThres[i] if tmpAr[j] > wigThres[i] else tmpAr[j] ) ))
			# end for j
		# end for chrm
	# end for i
	outFile.close()

def parseInputs( argv ):
	# Usage: python3.4 rna_chrm_maps.py [-ticks] [-o=out_id] [-c=chrms] [-b=bin_size] [-t=threshold] <gff_file> <fpkm_file> <wig_file> [wig_file]*

	outID = 'out'
	chrmList = None
	binSize = BINSIZE
	threshold = THRESHOLD
	percentile = PERCENTILE
	useTicks = False
	includeFileStr = None
	excludeFileStr = None
	startInd = 0
	
	for i in range(min(8,len(argv))):
		if argv[i].startswith('-o='):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-c='):
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith('-b='):
			try:
				binSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: bin size must be integer' )
				exit()
		elif argv[i].startswith('-t='):
			try:
				threshold = float(argv[i][3:])
				startInd += 1
			except ValueError:
				print( 'ERROR: threshold must be numeric' )
				exit()
		elif argv[i].startswith( '-p=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		elif argv[i] == '-ticks' or argv[i] == '-tick':
			useTicks = True
			startInd += 1
		elif argv[i].startswith( '-i=' ):
			if excludeFileStr != None:
				print( 'ERROR: cannot set include file and exclude file')
				exit()
			includeFileStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-e=' ):
			if includeFileStr != None:
				print( 'ERROR: cannot set include file and exclude file')
				exit()
			excludeFileStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a recognized option'.format( argv[i] ) )
			exit()
			
	
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedGraphFileAr = []
	for j in range(startInd+2,len(argv)):
		bedGraphFileAr+= [ argv[j] ]
	processInputs( gffFileStr, fpkmFileStr, bedGraphFileAr, outID, chrmList, binSize, threshold, percentile, useTicks, includeFileStr, excludeFileStr )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 rna_chrm_maps.py [-ticks] [-o=out_id] [-c=chrms] [-b=bin_size] [-t=threshold] [-i=include_file] [-e=exclude_file] <gff_file> <fpkm_file> [wig_file]*")
	else:
		parseInputs( sys.argv[1:] )
