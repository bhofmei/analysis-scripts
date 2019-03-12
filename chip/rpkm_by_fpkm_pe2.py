import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 rpkm_by_fpkm_pe.py [-strand] [-o=output_prefix] [-p=num_proc] [-l=labels] <gff_file> <fpkm_file> <bed_file> [bed_file]*
# computes ChIP RPKM for each gene from each bed file
# orders the output based on gene expression (FPKM)

NUMPROC=1

def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, outPre, numProc, sampleLabels, isStrand ):
	fileSamples = getSampleNames( bedFileStrAr )
	if sampleLabels == None:
		sampleNamesAr = fileSamples
	else:
		sampleNamesAr = sampleLabels
	
	info = '#from_script: rpkm_by_fpkm_pe.py; stranded:{:s}; gff_file:{:s}; fpkm_file:{:s}; files:{:s}'.format( str(isStrand), os.path.basename( gffFileStr ), os.path.basename( fpkmFileStr ), ' '.join( fileSamples ) )
	# read FPKM
	print( 'Reading FPKM' )
	fpkmAr = getFPKM( fpkmFileStr )
	# read GFF
	print( 'Reading GFF' )
	gffDict = readGFF( gffFileStr )
	
	# analyze bed files
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing {:d} samples with {:d} processors'.format( len(bedFileStrAr), numProc ) )
	results = [ pool.apply_async( processFiles, args=(bedFileStrAr[i], gffDict, isStrand ) ) for i in range(len(bedFileStrAr)) ]
	outDictAr = [ p.get() for p in results ]
	
	outFileStr = outPre + '_rpkm_by_fpkm.tsv'
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outDictAr, fpkmAr, sampleNamesAr, info )

def getSampleNames( fileStrAr ):
	
	names = []
	for fileStr in fileStrAr:
		tmpAr = fileStr.split( ',' )
		for i in range(len(tmpAr)):
			tmpF = os.path.basename( tmpAr[i] )
			rInd = tmpF.rfind('.')
			tmpAr[i] = tmpF[:rInd]
		names += ['/'.join(tmpAr)]
	return names
	
def getFPKM( fpkmFileStr ):
	'''
		returns the array or genes ordered by fpkm (highest first)
	'''
	# check if exists
	checkFPKM( fpkmFileStr )
	
	fpkmFile = open( fpkmFileStr + ".txt", 'r' )
	fpkmAr = []
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		i = lineAr[0].rfind('.')
		if line.startswith( 'Potri')==False and i != -1:
			lineAr[0] = lineAr[0][:i]
		fpkmAr += [ (lineAr[0],float(lineAr[1]))  ]
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

def processFiles( bedFileStrAr, gffDict, isStrand ):
	# check multiple files
	bedFiles = bedFileStrAr.split( ',' )
	print( 'Files in sample: {:s}'.format( ', '.join(bedFiles)))
	n = len(bedFiles)
	if n == 1:
		bedDict = processFile( bedFiles[0], gffDict, isStrand )
	else:	# loop through and combine
		bedDictAr = []
		for i in range(n):
			bedDictAr += [ processFile( bedFiles[i], gffDict, isStrand ) ]
		# combine
		bedDict = combineCounts( bedDictAr )
	return bedDict
		
def processFile( bedFileStr, gffDict, isStrand ):
	outDict = {}
	print( 'Processing', bedFileStr )
	bedDict, readCounts = readBed( bedFileStr, isStrand )
	
	# loop through genes
	for gene in gffDict.keys():
		info = gffDict[gene]
		chrm = info[0]
		start = info[1]
		end = info[2]
		strand = info[3]
		# compute RPKM
		if isStrand:
			chrm += strand
		counts = countsPerRegion( bedDict, chrm, start, end )
		adjusted = float( counts ) / ( (end-start+1) * readCounts ) * math.pow(10,9)
		outDict[gene] = adjusted
	return outDict

def readBed( bedFileStr, isStrand):
	'''
		creates a dictionary with each scaffold and those dictionaries are a dictionary
		for the positions of coordinates with frequency count from the bed file
		{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
		return the dictionary
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	readCount = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		readCount += 1
		lineAr =line.rstrip().split('\t')
		try:
			curScaf = lineAr[0]
			if isStrand:
				curScaf += lineAr[5]
			#pos = int( lineAr[1] ) + 1
			start = int( lineAr[1] ) + 1
			end = int( lineAr[2] ) + 1
			pos = int( math.floor( (end - start + 1)/2.0 + start ) )
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
	return bedDict, readCount

def countsPerRegion( bedDict, chrm, start, end ):
	counts = 0
	for i in range(start, end+1):
		try:
			dictEntry = bedDict.get(chrm).get(i)
			if dictEntry != None:
				counts += dictEntry
		except AttributeError:
			pass
	# end for
	return counts

def combineCounts( bedDictAr ):
	outDict = {}
	n = len(bedDictAr)
	for gene in bedDictAr[0].keys():
		vals = [ bedDict.get( gene ) for bedDict in bedDictAr ]
		if None in vals:
			continue
		avg = sum( vals ) / float( n )
		outDict[gene] = avg
	return outDict


def writeOutput( outFileStr, outDictAr, fpkmAr, sampleNamesAr, info ):
	headerAr = [ 'sample', 'gene.num', 'gene', 'rna.fpkm', 'chip.rpkm']
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' + '\t'.join(headerAr) + '\n' )
	numSamples = len(outDictAr)
	numGenes = len( fpkmAr )
	# loop through samples
	for i in range( numSamples ):
		# loop through genes
		for j in range( numGenes ):
			info = outDictAr[i].get( fpkmAr[j][0] )
			if info == None:
				continue
			# write output
			outStr = '{:s}\t{:d}\t{:s}\t{:.2f}\t{:f}\n'.format( sampleNamesAr[i], numGenes-j, fpkmAr[j][0], fpkmAr[j][1], info )
			outFile.write( outStr )
		# end for j
	# end for i
	outFile.close()

def parseInputs( argv ):
	outPre = 'out'
	numProc = NUMPROC
	sampleLabels = None
	isStrand = False
	startInd = 0
	
	for i in range(min(4,len(argv)-3)):
		if argv[i].startswith('-o='):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-l=' ):
			tmp = argv[i][3:]
			sampleLabels = tmp.split(',')
			startInd += 1
		elif argv[i] == '-strand':
			isStrand = True
			startInd += 1
		elif argv[i].startswith('-'):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	if sampleLabels != None and len(sampleLabels) != len( bedFileStrAr ):
		print( 'WARNING: number of labels does not match number of BED files. Ignoring labels.' )
		sampleLabels = None
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, outPre, numProc, sampleLabels, isStrand )


if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 rpkm_by_fpkm_pe.py [-strand] [-o=output_prefix] [-p=num_proc] [-l=labels] <gff_file> <fpkm_file> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
