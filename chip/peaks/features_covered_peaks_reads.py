import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3 features_covered_peaks_reads.py [-d] [-m=faction_covered] [-w=min_fraction_covered] [-o=out_id] [-s=promoter_size] [-f=fasta_index] [-r=repeat_gff] <gene_gff <peak_file>

UPFRAC=2/3
LOWFRAC=0.05
PROMOTSIZE=1000

def processInputs( gffFileStr, peakBedStr, repeatGffFile, fastaIndexFile, outID, upperFractionCovered, lowerFractionCovered, promoterSize, isDebug ):
	
	info = '#from_script:features_covered_peaks_reads.py; gff_file:{:s}; peak_bed_file:{:s}'.format( os.path.basename( gffFileStr ), os.path.basename( peakBedStr ) )
	print( 'Reading GFF' )
	geneDict, repeatDict, chrmDict = readGFF( gffFileStr )
	
	if repeatGffFile != None:
		print( 'Reading Repeat GFF' )
		repeatDict = readRepeatGFF( repeatGffFile )
		info += '; repeat_gff:{:s}'.format( os.path.basename( repeatGffFile ) )
	
	if fastaIndexFile != None:
		print( 'Reading FASTA index' )
		chrmDict = readFastaIndex( fastaIndexFile )
		info += '; fasta_index:{:s}'.format( os.path.basename( fastaIndexFile ) )
	
	# chromosome feature labeling
	print( 'Preparing for Peak Analysis' )
	featDict = prepFeatChrm( chrmDict, geneDict, repeatDict, promoterSize )
	info += '; promoter_size: {:d}; upper_fraction_threshold: {:.4f}; lower_fraction_threshold: {:.4f}'.format( promoterSize, upperFractionCovered, lowerFractionCovered )
	
	#genome heuristics
	gPer, pPer, tPer, iPer = computeGenomeHeuristics( featDict )
	info += '; gene_per: {:.2}; promot_per: {:.2f}; TE_per: {:.2f}; inter_per: {:.2f}'.format( gPer*100, pPer*100, tPer*100, iPer*100)
	print( 'Genome Breakdown:\n-genes: {:.2f}%\n-promoters: {:.2f}%\n-transposons: {:.2f}%\nintergenic: {:.2f}%'.format(gPer*100, pPer*100, tPer*100, iPer*100) )
	
	if outID == None:
		tmpN = os.path.basename( peakBedStr )
		rInd = tmpN.rfind( '.' )
		outID = tmpN[:rInd]
		
	# dealing with debug peak file
	if isDebug:
		debugFileStr = outID + '_peaks.tsv'
	else:
		debugFileStr = None
	
	# read peak file
	print( 'Reading and Analyzing Peak File' )
	outDict = readPeakFile( peakBedStr, featDict, upperFractionCovered, lowerFractionCovered, debugFileStr )
	
	# write output
	outFileStr = outID + '_out.tsv'
	print( 'Writing Output to {:s}'.format(outFileStr) )
	writeOutput( outDict, outFileStr, info )
	print( 'Done' )
	
def readGFF( gffFileStr ):
	'''
		reads GFF file which includes gene annotation and may or may not include
		repeat annotation
		returns 3 dictionaries: genes, repeats, chrm lengths
	'''
	geneDict = {}
	repeatDict = {}
	chrmDict = {}
	currChrm = None
	currMax = 0
	
	gffFile = open( gffFileStr, 'r' )
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		strand = lineAr[6]
		
		# if not gene, repeat, or TE, ignore
		if lineAr[2] not in [ 'gene', 'transposable_element_gene', 'transposon_fragment', 'transposable_element', 'similarity' ]:
			continue
		# new chromosome
		if currChrm != chrm:
			chrmDict = addToChrmDict( chrmDict, currChrm, currMax )
			currChrm = chrm
			currMax = 0
		# check max
		currMax = max( end, currMax )
		# add feature
		if lineAr[2] == 'gene':
			geneDict = addToDict( geneDict, chrm, (start, end, strand) )
		else:	# te, repeat
			repeatDict = addToDict( repeatDict, chrm, (start, end) )
	# end for line
	# add last chrm
	chrmDict = addToChrmDict( chrmDict, currChrm, currMax )
	gffFile.close()
	return geneDict, repeatDict, chrmDict

def addToChrmDict( chrmDict, chrm, inMax ):
	'''
		utility function that updates chrm length dictionary
		returns input dict with updated length if appropriate
	'''
	if chrm != None:
		prevCurrMax = chrmDict.get( chrm )
		if prevCurrMax == None:	# not previously known
			chrmDict[ chrm ] = inMax
		elif prevCurrMax < inMax:	# new value is larger
			chrmDict[ chrm ] = inMax
	return chrmDict
				
def addToDict( inDict, key, tup ):
	'''
		utility function that adds a tuple to dictionary
		returns input dict with tuple added to array for the key
	'''
	if inDict.get( key ) == None:
		inDict[key] = []
	inDict[key] += [ tup ]
	return inDict

def readRepeatGFF( repeatGffFile ):
	'''
		read GFF file that contains repeat annotations
		return dictionary; key is chrm name; value is array of tuples with
		(start, end) of repeats
	'''
	repeatDict = {}
	repeatFile = open( repeatGffFile, 'r' )
	
	for line in repeatFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		repeatDict = addToDict( repeatDict, chrm, (start, end) )
	# end for line
	repeatFile.close()
	return repeatDict

def readFastaIndex( fastaIndexStr ):
	'''
		read fasta index for chromosome lengths and names
		return dictionary; key is chrm name, value is length
	'''
	chrmDict = {}
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) length ...
		chrmDict[ lineAr[0] ] = int( lineAr[1] )
	fastaIndex.close()
	return chrmDict

def prepFeatChrm( chrmDict, geneDict, repeatDict, promoterSize ):
	'''
		takes in the annotation dictionaries, chrm lengths, and promoter sizes
		returns dictionary; key is chrm name, value is array of length 
		chrm-length
		each index in array is annotated "g", "i", "t", or "p" for gene, 
		intergenic, TE, and promoter respectively
	'''
	outDict = {}
	
	# loop through chrms
	for chrm in chrmDict.keys():
		tmpAr = [''] + ['i'] * chrmDict[chrm]	# assume intergenic
		# loop through repeats -> lower priority than genes
		if repeatDict.get( chrm ) != None:
			for repeat in repeatDict[chrm]:
				start = repeat[0]
				end = repeat[1]
				for i in range(start, end+1):
					tmpAr[i] = 't'
			# end loop through repeats
			
		# loop through genes -> higher priority
		# promoters have lower priority than genes, higher than repeats
		if geneDict.get( chrm ) == None:
			continue
		for gene in geneDict[chrm]:
			start = gene[0]
			end = gene[1]
			strand = gene[2]
			# promoter
			if strand == '+':
				pStart = max(1, start-promoterSize)
				for i in range(pStart, start):
					# overwrite promoter with intergenic/repeat space
					tmpAr[i] = ( 'p' if tmpAr[i] != 'g' else 'g' )
			else:
				pEnd = min(len(tmpAr)-1, end+promoterSize)
				for i in range(end+1, pEnd+1):
					tmpAr[i] = 'p'
			# gene
			for i in range(start, end+1):
				tmpAr[i] = 'g'
		# end loop through genes
		outDict[chrm] = tmpAr
	# end loop through chrm
	return outDict

def computeGenomeHeuristics( featDict ):
	'''
		compute proportion of genes, promoters, repeats/TEs, intergenic for the
		genome to get an idea of distribution
	'''
	gTotal = 0
	pTotal = 0
	tTotal = 0
	iTotal = 0
	total = 0
	
	for chrm in featDict.keys():
		featAr = featDict[chrm]
		gTotal += sum( [x == 'g' for x in featAr] )
		pTotal += sum( [x == 'p' for x in featAr] )
		tTotal += sum( [x == 't' for x in featAr] )
		iTotal += sum( [x == 'i' for x in featAr] )
		total += len( featAr ) - 1	# 0th index is empty
	perAr = [ float(x) / float(total) for x in [gTotal, pTotal, tTotal, iTotal] ]
	return perAr[0], perAr[1], perAr[2], perAr[3]

def readPeakFile( peakBedStr, featDict, upperFractionCovered, lowerFractionCovered, peaksDebugStr ):
	outCountDict = { 'g':0,'p':0,'t':0,'i':0,'o':0,'a':0}
	chrms = list( featDict.keys() )
	
	if peaksDebugStr != None:
		peaksOutFile = open( peaksDebugStr, 'w' )
		headerAr = ['chrm','start','end','name','assignment','g','i','p','t']
		peaksOutFile.write( '#{:s}\n'.format( '\t'.join( headerAr ) ) )
	
	peakFile = open( peakBedStr, 'r' )
	for line in peakFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		chrm = lineAr[0]
		if chrm not in chrms:
			continue
		start = int( lineAr[1] ) + 1
		end = int( lineAr[2] ) + 1
		peakDict = {'g':0,'p':0,'t':0,'i':0}
		# loop through peak
		for i in range(start, end+1):
			val = featDict[chrm][i]
			peakDict[ val ] += 1
		peakType = definePeakType( peakDict, (end-start+1), upperFractionCovered, lowerFractionCovered)
		
		# debuging if necessary
		if peaksDebugStr != None:
			tmpVals = [ str(peakDict[x]) for x in sorted(peakDict.keys()) ]
			outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:s}\n'.format( chrm, start, end, lineAr[3], peakType, '\t'.join( tmpVals ) )
			peaksOutFile.write( outStr )
		# update counter
		outCountDict[peakType] += 1
		outCountDict['a'] += 1
		
	peakFile.close()
	if peaksDebugStr != None:
		peaksOutFile.close()
	return outCountDict
			
def definePeakType( inDict, peakLen, upFrac, lowFrac ):
	typesAr = ['g', 'p', 't', 'i']
	valsAr = [ float( inDict[typ] ) / peakLen for typ in typesAr ]
	
	# no features above upper threshold -> overlap/ambiguous
	if sum ( [ x >= upFrac for x in valsAr ] ) == 0:
		return 'o'
	# more than one feature above low threshold -> overlap
	elif sum ( [ x > lowFrac for x in valsAr ] ) > 1:
		return 'o'
	# else it has a type -> assign it that type
	else:
		mInd = valsAr.index( max(valsAr) )
		return typesAr[mInd]

def writeOutput( outDict, outFileStr, info ):
	'''
		function to write output to the file outFileStr
	'''
	outFile = open( outFileStr, 'w' )
	headerAr = [ 'peakType','count','percent' ]
	outFile.write( info + '\n' + '\t'.join(headerAr) + '\n')
	
	totalPeaks = outDict['a']
	peakTypes = ['g','p','t','i','o','a']
	peakTypeFormat = [ 'genes', 'promoters', 'transposons/repeats', 'intergenic', 'overlap', 'all' ]
	
	for i in range(len(peakTypes)):
		val = outDict[ peakTypes[i] ]
		outStr = '{:s}\t{:d}\t{:.2f}\n'.format( peakTypeFormat[i], val, float(val) / totalPeaks * 100 )
		outFile.write( outStr )
	# end for

def parseInputs( argv ):
	upperFractionCovered = UPFRAC
	lowerFractionCovered = LOWFRAC
	outID = None
	promoterSize = PROMOTSIZE
	fastaIndexFile = None
	repeatGffFile = None
	isDebug = False
	startInd = 0
	
	for i in range(min(7,len(argv))):
		if argv[i].startswith( '-m=' ):
			try:
				upperFractionCovered = float( argv[i][3:] )
				if upperFractionCovered > 1:
					upperFractionCovered = upperFractionCovered / 100
				startInd += 1
			except ValueError:
				print( 'ERROR: fraction covered parameter must be numeric' )
				exit()
		elif argv[i].startswith( '-w=' ):
			try:
				lowerFractionCovered = float( argv[i][3:] )
				if lowerFractionCovered > 1:
					lowerFractionCovered = lowerFractionCovered / 100
				startInd += 1
			except ValueError:
				print( 'ERROR: lower fraction covered parameter must be numeric' )
				exit()
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			fastaIndexFile = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-r=' ):
			repeatGffFile = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-s=' ):
			try:
				promoterStr = argv[i][3:]
				if promoterStr.endswith( 'k' ) or promoterStr.endswith( 'K' ):
					promoterSize = int( promoterStr[:-1] ) * 1000
				elif promoterStr.endswith( 'm' ) or promoterStr.endswith( 'M' ):
					promoterSize = int( promoterStr[:-1] ) * 1000000
				else:
					promoterSize = int( promoterStr )
				startInd += 1
			except ValueError:
				print( 'ERROR: promoter size must be an integer' )
				exit()
		elif argv[i] == '-d':
			isDebug = True
			startInd += 1
		elif argv[i] in [ '-h', '-help', '--help' ]:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for argv
	gffFileStr = argv[startInd]
	peakBedStr = argv[startInd+1]
	
	processInputs( gffFileStr, peakBedStr, repeatGffFile, fastaIndexFile, outID, upperFractionCovered, lowerFractionCovered, promoterSize, isDebug )

def printHelp():
	print( 'Usage: python3 features_covered_peaks.py [-d] [-m=faction_covered] [-w=min_fraction_covered] [-o=out_id] [-s=promoter_size] [-f=fasta_index] [-r=repeat_gff] <gene_gff <peak_file>' )
	print( 'Determines which peaks cover different types of features;\nfeatures include: genes, promoters, repeats/TEs, intergenic;\npeaks that overlap multiple types are ignored' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
