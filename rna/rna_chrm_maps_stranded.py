import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 rna_chrm_maps_stranded.py [-ticks] [-o=out_id] [-b=bin_size] [-t=threshold] [-m=min_fpkm] [-p=percentile] [-c=wt_pos] <gff_file> <fpkm_file> [bedgraph_file]*
# produce information to create genome-wide maps of genes and strand-specific
# up/down regulation of genes/anti-sense RNA

BINSIZE=100
THRESHOLD=1
TICKVAL=0
MINFPKM=0
PERCENTILE=0.95

def processInputs( gffFileStr, fpkmFileStr, bedFileAr, outID, chrmList, binSize, threshold, minFPKM, useTicks, basePos, percentile ):
	info = "#from_script:rna_chrm_maps_stranded.py; bin_size:{:d}; threshold:{:.1f}; min_fpkm:{:.1f}; use_ticks:{:s}; gff_file:{:s}; fpkm_file:{:s}".format( binSize, threshold, minFPKM, str(useTicks), os.path.basename(gffFileStr), os.path.basename(fpkmFileStr) )
	# read gff file -> get strand and coordinates
	print( 'reading gff...' )
	chrmDict, gffDict, geneDict = readGFF( gffFileStr )
	if chrmList == None:
		chrmList = list(chrmDict.keys())
	print(len(geneDict.keys()))
	#print( chrmDict )
	# process genes
	print('processing genes')
	geneValDict= processGenes( gffDict, chrmDict, chrmList, binSize )
	# write file and delete dictionaries
	outGeneStr = outID + '_genes.tsv'
	print('writing gene information to {:s}...'.format(outGeneStr))
	writeGenes(outGeneStr, geneValDict, info)
	del gffDict, geneValDict
	
	# read fpkm
	print('reading fpkm file...')
	sampleNameAr, fpkmDict = readFPKM( fpkmFileStr, list(geneDict.keys()) )
	print( 'base sample:', sampleNameAr[ basePos ] )
	print('processing fpkm...' )
	fpkmValDict = processFPKM( fpkmDict, chrmDict, chrmList, binSize, sampleNameAr, geneDict, minFPKM, basePos )
	print( 'cleaning fpkm...' ) 
	if useTicks:
		fpkmCleanDict = cleanFPKMTicks( fpkmValDict, threshold )
	else:
		fpkmCleanDict = cleanFPKM( fpkmValDict, threshold )
	outFPKMStr = outID + '_fpkm'+('_ticks' if useTicks else '') +  '.tsv'
	print('writing fpkm information to {:s}...'.format(outFPKMStr))
	writeFPKM( outFPKMStr, fpkmCleanDict, sampleNameAr, info, useTicks )
	
	# handle bedgraph files
	if len(bedFileAr) > 0:
		bedNames = [ getSampleName(x) for x in bedFileAr ]
		bedFiles = [os.path.basename( x ) for x in bedFileAr ]
		infoBed = info + ';bedgraph_files:{:s};percentile:{:.3f}'.format( ','.join( bedFiles), percentile )
		print('reading bedgraph files...')
		bedDictAr = [ readBedGraph( x, chrmDict, chrmList, binSize ) for x in bedFileAr ]
		bedThres = [ correctBed(w, percentile) for w in bedDictAr ]
		outBedStr =  outID + '_{:d}_bedgraph.tsv'.format(int(percentile*100))
		print('writing bedgraph information to {:s}...'.format(outBedStr))
		writeBed( outBedStr, bedDictAr, bedNames, bedThres, infoBed )
	
def readGFF( gffFileStr ):
	gffFile = open( gffFileStr, 'r' )
	chrmDict = {}
	geneDict = {}
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
		# for gene maps
		elif lineAr[2] == "gene":
			chrm = lineAr[0]
			start = int(lineAr[3])
			end = int(lineAr[4])
			#name = getGeneName( lineAr[8] )
			if gffDict.get(chrm) == None:
				gffDict[chrm] = []
			gffDict[chrm] += [( start, end,lineAr[6] )]
			# in case chromosomes aren't part of GFF
			if chrm not in chrmSet:
				chrmSet.add(chrm)
				chrmDictGene[chrm] = end
			if end > chrmDictGene[chrm]:
				chrmDictGene[chrm] = end
		# for fpkms
		elif lineAr[2] == "mRNA":
			rName = getRNAName( lineAr[8] )
			start = int(lineAr[3])
			end = int(lineAr[4])
			strand = lineAr[6]
			geneDict[rName] = (chrm, start, end, strand)
	# end for
	if chrmDict == {}:
		chrmDict = chrmDictGene
	gffFile.close()
	return chrmDict, gffDict, geneDict

def getRNAName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find('-')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def processGenes( gffDict, chrmDict, chrmList, binSize ):
	outDict = {}
	# loop through chromosomes
	for chrm in chrmList:
		numBins = int(math.ceil(chrmDict[chrm]/float(binSize)))
		outDict[chrm] = [0] * numBins
		# loop through genes
		geneAr = gffDict.get(chrm)
		if geneAr == None:
			continue
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

def writeGenes( outFileStr, outDict, info ):
	outFile = open( outFileStr, 'w' )
	header = info + '\n' + "#chrm\tpos\tvalue\tstrand\n"
	outFile.write(header)
	# loop through chromosomes
	for chrm in sorted(outDict.keys()):
		# loop through bins
		binAr = outDict[chrm]
		strandAr = ["NA"] * len( binAr )
		for j in range(len(binAr)):
			if binAr[j] >= 1:
				binAr[j] = 1
				strandAr[j] = 'pos'
			elif binAr[j] <= -1:
				binAr[j] = -1
				strandAr[j] = 'neg'
			outFile.write( '{:s}\t{:d}\t{:d}\t{:s}\n'.format( chrm, j, binAr[j], strandAr[j] ) )
	# end for chrm
	outFile.close()
	
def readFPKM( fpkmFileStr, geneAr ):
	sampleNameAr = []
	indexAr = []
	classIndex = 1
	geneIndex = 4
	isHeader = True
	fpkmDict = {}
	
	fpkmFile = open( fpkmFileStr, 'r' )
	for line in fpkmFile:
		lineAr = line.rstrip().split('\t')
		# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id
		# (4) gene_short_name (5) gene_description (6) tss_id (7) locus 
		# (8) length (9) coverage (10) sName_FPKM (11) sName_conf_lo 
		#(12) sName_conf_hi (13) sName_status
		if isHeader:
			# header
			for j in range(len(lineAr)):
				if lineAr[j] == "gene_short_name":
					geneIndex = j
				elif lineAr[j] == "class_code":
					classIndex = j
				elif lineAr[j].endswith('_FPKM'):
					indexAr += [j]
					sampleNameAr += [ lineAr[j].replace('_FPKM','')]
			isHeader = False
		else:
			geneName = lineAr[geneIndex]
			if geneName not in geneAr:
				continue
			# gene name need to include class
			if lineAr[classIndex] == "x":
				locusName = "~" + geneName
			else:
				locusName = geneName
			
			fpkmValAr = [0]*len(indexAr)
			for k in range(len(indexAr)):
				if lineAr[indexAr[k]+3]=="OK":
					fpkmValAr[k] +=  float( lineAr[indexAr[k]] ) 
			# add to dictionary
			fpkmDict[locusName] = fpkmValAr
	# end for
	fpkmFile.close()
	return sampleNameAr, fpkmDict
	
def processFPKM( fpkmDict, chrmDict, chrmList, binSize, sampleNamesAr, geneDict, minFPKM, basePos ):
	fpkmValDict = {}
	for chrm in chrmList:
		numBins = int(math.ceil(chrmDict[chrm]/float(binSize)))
		# two strands
		chrm1 = chrm + "-"	# negative strand
		chrm2 = chrm + "+"  # positive strand
		fpkmValDict[chrm1] = [ [0]*numBins for i in range(len(sampleNamesAr)) ]
		fpkmValDict[chrm2] = [ [0]*numBins for i in range(len(sampleNamesAr)) ]
	
	# loop through genes
	for locus in fpkmDict.keys():
		# get fpkm values
		fpkmAr = fpkmDict[locus]
		# check anti-sense
		if locus.startswith("~"):
			isAnti = True
			locus = locus[1:]
		else:
			isAnti = False
		# get gene information
		tup = geneDict.get( locus )
		if tup == None:
			continue
		#print(locus, tup )
		chrm, start, end, strand = tup
		# get bin information
		bStart = start // binSize
		bEnd = end // binSize
		# check fpkm values -> at least one greater than min fpkm
		t = [ (1 if x < minFPKM else 0) for x in fpkmAr ]
		ts = sum(t)
		if ts == len(fpkmAr):
			adj = 0
		else:
			adj = 1
		vChrm = getAdjustChrm( chrm, strand, isAnti )
		#print( "-", ts, vChrm )
		# loop through sample comparison
		for j in range(len(fpkmAr)):
			val = computeChange(fpkmAr[basePos],fpkmAr[j])
			for bin in range(bStart, bEnd + 1):
				fpkmValDict[vChrm][j][bin] += adj * val
		
	# end for genes
	return fpkmValDict

def getAdjustChrm( chrm, strand, isAnti ):
	if isAnti == False:
		return chrm + strand
	elif strand == "-":
		return chrm +  "+"
	else:
		return chrm + "-"

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

def cleanFPKMTicks( outDict, threshold ):
	newDict = {}
	# loop through chms
	for chrm in sorted(outDict.keys()):
		newDict[chrm] = [ None ] *len(outDict[chrm])
		# loop through samples
		for i in range(len(outDict[chrm])):
			arr = outDict[chrm][i]
			newAr = ["NA"] * len(arr)
			# loop through bins:
			for j in range(len(arr)):
				# check upper threshold
				if arr[j] >= threshold:
					#print( j, max(0,j-2), min(len(arr),j+3) )
					for k in range( max(0,j-TICKVAL), min(len(arr), j+TICKVAL+1)):
						newAr[k] = "up"
				# check lower threshold
				elif arr[j] <= -1 * threshold:
					#print( j, max(0,j-2), min(len(arr),j+3) )
					for k in range( max(0,j-TICKVAL), min(len(arr), j+TICKVAL+1)):
						newAr[k] = "down"
						
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

def writeFPKM( outFileStr, outDict, sampleAr, info, useTicks ):
	outFile = open( outFileStr, 'w' )
	header = info + '\n'+'#chrm\tsample\tpos\tvalue\tstrand\n'
	outFile.write(header)
	
	# loop through chromosomes
	for chrm in sorted(outDict.keys()):
		if chrm.endswith( '-' ):
			strand = 'neg'
		else:
			strand = 'pos'
		cChrm = chrm[:-1]
		# loop through samples
		for i in range(len(sampleAr)):
			tmpAr = outDict[chrm][i]
			#print(tmpAr)
			# loop through bins
			for j in range(len(tmpAr)):
				if useTicks:
					#print(cChrm, sampleAr[i], j, tmpAr[j], strand)
					outFile.write('{:s}\t{:s}\t{:d}\t{:s}\t{:s}\n'.format( cChrm, sampleAr[i], j, tmpAr[j], strand ))
				else:
					outFile.write('{:s}\t{:s}\t{:d}\t{:.4f}\t{:s}\n'.format( cChrm, sampleAr[i], j, tmpAr[j], strand ))
			# end for j
		# end for i
	# end for chrm
	outFile.close()
	
def getSampleName( fileStr ):
	# bed file
	leftIndex = fileStr.rfind('/')
	rightIndex = fileStr.rfind('.')
	sampleName = fileStr[leftIndex+1:rightIndex]
	return sampleName

def readBedGraph( fileStr, chrmDict, chrmList, binSize ):
	bedDict = {}
	bedFile = open( fileStr, 'r' )

	# set up bedDict
	for chrm in chrmList:
		n = int(math.ceil(chrmDict[chrm]/float(binSize)))
		bedDict[chrm] = [0]*n

	for line in bedFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) start (2) end (3) value
		# first time seeing the chrm
		chrm = lineAr[0]
		start = int( lineAr[1] ) + 1
		#end = int( lineAr[2] ) + 1
		value = float( lineAr[3] )
		bin = int( start ) // binSize
		if bedDict.get( chrm ) != None:
			bedDict[chrm][bin] += value
	bedFile.close()
	return bedDict

def correctBed( bedDict, percentile ):
	numList = []
	# loop through chromosomes and collect numbers
	for chrm in bedDict.keys():
		numList += bedDict[chrm]
	if percentile == 1:
		return max(numList)
	numList.sort()
	ind = math.ceil(percentile * len( numList ) - 1)
	try:
		return numList[ind]
	except IndexError:
		return numList[-1]

def writeBed( outFileStr, bedDictAr, bedNameAr, bedThres, info ):
	outFile = open( outFileStr, 'w' )
	header = info + '\n' + '#bed\tchrm\tpos\tvalue\n'
	outFile.write(header)
	#print(wigDictAr)
	# loop through samples
	for i in range(len(bedDictAr)):
		# loop through chrms
		for chrm in sorted(bedDictAr[i].keys()):
			tmpAr = bedDictAr[i][chrm]
			# loop through bins
			for j in range(len(tmpAr)):
				outFile.write('{:s}\t{:s}\t{:d}\t{:.4f}\n'.format( bedNameAr[i], chrm, j, ( bedThres[i] if tmpAr[j] > bedThres[i] else tmpAr[j] ) ))
			# end for j
		# end for chrm
	# end for i
	outFile.close()

def parseInputs( argv ):
	binSize = BINSIZE
	threshold = THRESHOLD
	chrmList = None
	minFPKM = MINFPKM
	outId = 'out'
	percentile=PERCENTILE
	useTicks = False
	basePos = 0
	startInd = 0
	
	for i in range(min(8,len(argv))):
		if argv[i].startswith('-o='):
			outId = argv[i][3:]
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
		elif argv[i].startswith( '-m=' ):
			try:
				minFPKM = float( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: min fpkm must be numeric' )
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
		elif argv[i] == '-ticks':
			useTicks = True
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			try:
				basePos = int( argv[i][4:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: wild type position must be integer' )
				exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
			
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	wigFileAr = []
	
	for j in range(startInd+2,len(argv)):
		wigFileAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, wigFileAr, outId, chrmList, binSize, threshold, minFPKM, useTicks, basePos, percentile )

def printHelp( ):
	print ("Usage: python3 rna_chrm_maps_stranded.py [-ticks] [-o=out_id] [-b=bin_size] [-c=chrm_list] [-t=threshold] [-m=min_fpkm] [-c=wt_pos] [-p=percentile] <gff_file> <fpkm_file> [wig_file]*")
	print( 'Produces data for RNA-seq fold-change chromosome plots;\nplots show plus and minus strand separately' )
	print( 'Required: ')
	print( 'gff_file\tpath to GFF format annotation file of genes' )
	print( 'fpkm_file\tpath to cuffdiff output of sample FPKM values' )
	print( 'wig_file\tpath to wig file' )
	print( 'Optional:' )
	print( '-ticks\t\tuse tick marks not fold-change dependent coloring' )
	print( '-o=out_id\tstring identifier for output file [default: out]' )
	print( '-b=bin_size\tbin size to break chromosome into [default 100]' )
	print( '-c=chrm_list\tcomma-separated list of chromosomes to include\n\t\t[default all in GFF]' )
	print( '-t=log2_threshold\twith ticks, genes above this are marked\n\t\twithout ticks, max fold-change value to improve vizualization\n\t\t[default 1]' )
	print( '-m=min_fpkm\tminimum fpkm value; at least one sample must be above this\n\t\t[default 0]' )
	print('-p=percentile\tpercentile correction for wig files [default 0.95]' )
	print( '-c=wt_pos\t0-based index indicating the wild-type/control sample when\n\t\tcomputing fold-changes [default 0]' )
	
if __name__ == "__main__":
	if len(sys.argv) < 3:
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
