import sys, math, glob, multiprocessing, subprocess, os, bisect
from bioFiles import *

# Usage: python3.4 heatmap_from_ortho_bed_allc_fr_pe.py [-1[,2] | -2[,1]] [-a|-r|-g|-n] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] [-f=ortholog_file] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*
# produce heatmaps of BED/AllC files based on ortholog options 
# note: gbm status always decided by the gbm status of the ortholog listed first

NUMPROC=2
NUMBINS=20
STREAMSIZE=1000
PERCENTILE=.95

ORTHO_FILE='/Users/bhofmeister/Documents/Research/collaborations/cmt3/arabidopsis_compare/Esalsugineum_Athaliana_bin.out'
#ORTHO_FILE='/escratch4/bth29393/bth29393_Jun_02/Ath_Esa_pairs.txt'

def processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, genesIncluded, orthoFile, speciesInt, fpkmInt ):
	''' 
	'''
	sampleNamesAr = getSampleNames( bedFileStrAr )
	print( 'Number of bins per gene: {:d}\nNumber of bins upstream/downstream: {:d}\nUpstream/Downstream size: {:d}\nGenes included: {:s}\nSpecies: {:d}\nFPKM Species: {:d}\nPercentile for correction: {:.3f}\nSamples included: {:s}\nGGF file: {:s}\nFPKM file: {:s}\n'.format( numBins, numBinsStream, streamSize, genesIncluded, speciesInt, fpkmInt, percentile, ', '.join(sampleNamesAr), os.path.basename(gffFileStr) , os.path.basename(fpkmFileStr) ) )
	chipBedStr = [ os.path.basename( x ) for x in bedFileStrAr ]
	info = '#from_script:heatmap_from_ortho_bed_allc_fr_pe.py; numBins:{:d}; numBinsStream:{:d}; streamSize:{:d}; percentileCorrection:{:.3f}; ortho_file:{:s}; gff_file:{:s}; fpkm_file:{:s}; chip_bed_files:{:s}'.format( numBins, numBinsStream, streamSize, percentile, os.path.basename(orthoFile), os.path.basename(gffFileStr) , os.path.basename(fpkmFileStr), ','.join(chipBedStr) )
	
	if genesIncluded == 'all':
		subsetD = None
	else:
		print( 'Reading ortholog file' )
		subsetD = readGBMPairs( genesIncluded, orthoFile, speciesInt, fpkmInt )
		# subsetDict is dictionary with fpkm gene as key and report gene as value
	print( 'Reading FPKM file' )
	#print( list(subsetD.keys())[0:3] )
	fpkmFile = FileFPKM( fpkmFileStr )
	fpkmAr, fpkmValAr = fpkmFile.getFPKMValueArray( subsetDict=subsetD )
	print( 'LEN:', len(fpkmAr), ',', len(fpkmValAr) )
	
	print( 'Reading GFF file' )
	gffFile = FileGFF( gffFileStr )
	geneDict, chrmFormat = gffFile.getGeneDict( chrm=True, rmPeriod=True )
	#print( list(geneDict.keys())[0:5] )
	if genesIncluded == 'all':
		outFileStrAr = [ '{:s}_{:d}_{:d}_{:d}.tsv'.format( sampleName, streamSize, numBins, int(percentile*100) ) for sampleName in sampleNamesAr ]
	else:
		outFileStrAr = [ '{:s}_{:d}_{:d}_{:d}_{:s}.tsv'.format( sampleName, streamSize, numBins, int(percentile*100), genesIncluded ) for sampleName in sampleNamesAr ]
	pool = multiprocessing.Pool( processes=numProc )
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	results = [ pool.apply_async( processFile, args=(bedFileStrAr[i], outFileStrAr[i], fpkmAr, geneDict, numBins, numBinsStream, streamSize, percentile, chrmFormat, fpkmValAr, info )) for i in range( len(bedFileStrAr) ) ]
	
	fin = [ p.get() for p in results ]
	if sum(fin) != len( bedFileStrAr ):
		print( 'ERROR: not enough files written' )
		exit()
	print( 'Done' )

def readGBMPairs( genesIncluded, orthoFile, speciesInt, fpkmInt ):
	
	outDict = {}
	inFile = open( orthoFile, 'r' )
	ind = speciesInt * 2 + -2
	fInd = fpkmInt * 2 + -2
	for line in inFile:
		lineAr = line.rstrip().split()
		i = lineAr[0].rfind('.')
		if i != -1:
			lineAr[0] = lineAr[0][:i]
		# (0) athaliana cds (1) UM/CG (2) eutrema cds (3) UM/CG
		try:
			# ortholog
			if genesIncluded == 'orthologs':
				outDict[ lineAr[fInd] ] = lineAr[ind] 
			# gbm
			elif genesIncluded == 'gBM_orthologs' and lineAr[1] == 'CG':
				outDict[ lineAr[fInd] ] = lineAr[ind] 
			# non-gBM
			elif genesIncluded == 'non_gBM_orthologs' and lineAr[1] == 'UM':
				outDict[ lineAr[fInd] ] = lineAr[ind] 
		except IndexError:
			print( 'ERROR:', line )
	inFile.close()
	return outDict

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

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

def processFile( fileStr, outFileStr, fpkmAr, gffDict, numBins, numBinsStream, streamSize, percentile, chrmForm, fpkmValAr, info ):
	isBed = checkFile( fileStr )
	
	if isBed:
		bedFile = FileBED_FR( fileStr )
		print( 'Processing {:s}'.format( str(bedFile) ) )
		bedDict, counts = bedFile.getBedDict()
		outMatrix, numList = processBed( bedDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize, counts )
		thresh = calculatePercentile( percentile, numList )
		print( 'Writing output for {:s}'.format( outFileStr ) )
		writeOutput( outFileStr, outMatrix, thresh, fpkmAr, fpkmValAr, info )
		return 1
	
	else:
		print( 'Processing {:s}'.format( fileStr ) )
		# this file has all chromosomes and we only want CG methylation
		allcFile = FileAllC_full( fileStr )
		allCDictAll = allcFile.getAllCDict( chrmFormat=chrmForm)
		allCDict = allCDictAll['CG']
		
		outMatrix, numList = processAllC( allCDict, fpkmAr, gffDict, numBins, numBinsStream, streamSize )
		#thresh = calculatePercentile( percentile, numList )
		thresh=1
		print( 'Writing output for {:s}'.format( outFileStr ) )
		writeOutput( outFileStr, outMatrix, thresh, fpkmAr, fpkmValAr, info )
		return 1

def checkFile( fileStr ):
	if fileStr.endswith( '.bed' ):
		return True
	else:
		# allC file
		return False

def processBed( bedDict, fpkmAr, geneDict, numBins, numBinsStream, streamSize, bpCount ):
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
			print( 'ERROR: no gene', gene )
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
		outAr = [ x / bpCount * 1000000 for x in outAr ]
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
			outAr[i] = (float(methAr[i][0]) / float(methAr[i][1]))
	return outAr

def writeOutput( outFileStr, outMatrix, thresh, fpkmAr, fpkmValAr, info ):

	outFile = open( outFileStr, 'w' )
	#geneCount = len( outMatrix )
	headerAr = [ 'geneNum','geneName', 'fpkm', 'bin', 'value' ]
	outFile.write( info + '\n' + '\t'.join( headerAr ) + '\n' )
	n=len(outMatrix)
	for i in range( len(outMatrix) ):
		vals = outMatrix[i]
		if vals[0] == -1:
			print( 'skipping: {:s}, no information'.format( fpkmAr[i] ) )
			continue
		for j in range(len(vals)):
			if vals[j] > thresh:
				vals[j] = thresh
			outStr = '{:d}\t{:s}\t{:s}\t{:d}\t{:f}\n'.format( n-i, fpkmAr[i], fpkmValAr[i], j, vals[j] )
			outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	numProc = NUMPROC
	numBins = NUMBINS
	numBinsStream = NUMBINS
	streamSize = STREAMSIZE
	percentile = PERCENTILE
	genesIncluded = 'all'
	speciesInt = 0
	fpkmInt = 0
	startInd = 0
	orthoFile = ORTHO_FILE
	
	for i in range( min(12, len(argv) ) ):
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
		elif argv[i].startswith('-1'):
			if speciesInt != 0:
				print( 'ERROR: cannot specify -1 with -2' )
				exit
			speciesInt = 1
			if argv[i] == '-1,2':
				fpkmInt = 2
			startInd += 1
		elif argv[i].startswith('-2'):
			if speciesInt != 0:
				print( 'ERROR: cannot specify -2 with -1' )
				exit()
			speciesInt = 2
			if argv[i] == '-2,1':
				fpkmInt = 1
			startInd += 1
		elif argv[i].startswith( '-h' ):
			printHelp()
			exit()
		elif argv[i].startswith ( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	
	if genesIncluded == 'all' and speciesInt != 0:
		print( 'ERROR: do not specify species int when using all genes' )
		exit()
	elif speciesInt == 0 and genesIncluded != 'all':
		speciesInt = 2
	
	if fpkmInt == 0:
		fpkmInt = speciesInt
	
	gffFileStr = argv[startInd]
	fpkmFileStr = argv[startInd+1]
	bedFileStrAr = []
	
	for j in range( startInd + 2, len(argv) ):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( gffFileStr, fpkmFileStr, bedFileStrAr, numProc, numBins, numBinsStream, streamSize, percentile, genesIncluded, orthoFile, speciesInt, fpkmInt )
	
def printHelp():
	print('Usage: python3.4 heatmap_from_ortho_bed_allc_fr_pe.py [-1[,2]|-2[,1]] [-a|-r|-g|-n] [-p=num_proc] [-b=num_bins[,num_bins_stream]] [-s=stream_size] [-t=percentile] [-f=ortholog_file] <gff_file> <fpkm_file> <bed_file> [bed_file | allc_file]*\n')
	print('Required:\ngff_file\tgff file of species interested in')
	print('fpkm_values\tfpkm value file of species to order genes by' )
	print( 'bed_file\tBED file of ChIP interested in' )
	print( 'allc_file\tallC file of methylation; all chrms in one file\n')
	print( 'Optional:\n-1/2\tspecify species for genes/fpkm values as listed in ortholog file' )
	print( '\t\t-1  use genes and fpkm values listed first in ortholog file\n\t\t-1,2  use genes listed first in ortholog file but fpkm values are for genes listed second\n\t\t-2  use genes and fpkm values listed second in ortholog file\n\t\t-2,1  use genes listed second in ortholog file but fpkm values are for genes listed first' )
	print( '-a/r/g/n\tspecify which type of genes to include\n' )
	print( '\t\t-a  include all genes\n\t\t-r  include only orthologs\n\t\t-g  include only gBM orthologs\n\t\t-n  include only non-gBM orthologs' )
	print( '-f=ortholog_file\tuse this ortholog file instead' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )
	print( '-b=N[,W]\tuse N bins for gene body and W bins for upstream/downstream [default 20,20]' )
	print( '-s\tnumber of base pairs to include upstream and downstream [default 1000]' )
	print( '-t=D\tpercentile to correct largest values by' )
	print( '-h\tprint this help screen' )
	

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
