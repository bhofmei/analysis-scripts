import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
import pysam #, pybedtools
from functools import reduce

# Usage: python heatmap_abs_bam_pe.py [-h] [-q] [-z] [-n] [-b=num_bins] [-d=distance]
# [-o=out_id] [-t=percentile] [-r=gene_order] <gff_file> <bam_file> [bam_file]*

NUMBINS=10
DISTANCE=200
NUMPROC=1
PERCENTILE = 0.99

def processInputs( gffFileStr, bamFileAr, numProc, numBins, distance, outId, percentile, geneOrderStr, labels, isPrint, isCompress, isScaled ):
	if isPrint:
		print( 'GFF File:', os.path.basename(gffFileStr) )
		print( 'BAM files:', ', '.join([os.path.basename(file) for file in bamFileAr]) )
		print( 'Number of bins:', numBins )
		print( 'Distance:', distance )
		print( 'Percentile:', percentile )
		if geneOrderStr != None:
			print( 'Gene order file:', os.path.basename(geneOrderStr) )
	
	# get sample names
	if labels == None:
		labels = [ getSampleName(file) for file in bamFileAr ]
	
	# get gene array
	if isPrint:
		print( 'Getting genes from', os.path.basename(gffFileStr) )
	geneDict, geneLenAr = readGFF( gffFileStr )
	
	if geneOrderStr == None:
		# gene order by length
		geneLenAr.sort(reverse = True)
		geneOrder = [ g[1] for g in geneLenAr ]
		if isPrint:
			print( 'Ordering genes by length' )
	else:
		if isPrint:
			print( 'Ordering genes as defined in', os.path.basename(geneOrderStr) )
		geneOrder = readGeneOrder( geneOrderStr )
	
	
	info = 'from_script: heatmap_abs_bam_pe.py; gff_file: '
	info += os.path.basename(gffFileStr)
	if geneOrderStr != None:
		info += '; gene_order_file: ' + os.path.basename(geneOrderStr)
	info += '; distance: ' + str(distance) + '; num_bins: ' + str(numBins)
	info +=  '; is_scaled: ' + str(isScaled) + '; percentile: ' + str(percentile) 
	
	if isPrint:
		print( 'Begin processing {:d} files with {:d} processors'.format( len(bamFileAr), numProc) )
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(bamFileAr[i], labels[i], geneDict, geneOrder, distance, numBins, percentile, outId, isCompress, isScaled, isPrint, info) ) for i in range( len(bamFileAr) ) ]
	outN = [ p.get() for p in results ]
	print( 'Done' )

def getSampleName( fileStr ):
	bName = os.path.basename(fileStr)
	rIndex = bName.rfind('.')
	if rIndex != -1:
		return bName[:rIndex]
	else:
		return bName
		
def readGFF( gffFileStr ):
	
	gffFile = open( gffFileStr, 'r' )
	
	gffDict = {}
	geneLenAr = []
	
	for line in gffFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or len(lineAr) < 9:
			continue
		# only care about genes
		if lineAr[2] != 'gene':
			continue
		
		chrm = lineAr[0]
		start = int(lineAr[3])-1
		end = int(lineAr[4])
		strand = lineAr[6]
		gLength = end - start
		gName = getGeneName( lineAr[8] )
		gffDict[gName] = (chrm, start, end, strand)
		geneLenAr += [(gLength, gName)]
	# end for line
	
	gffFile.close()
	return gffDict, geneLenAr
	
def getGeneName( notesStr ):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	n = notesStr[adIndex:]
	if endIndex != -1:
		n = notesStr[adIndex:endIndex+adIndex]
	return n	

def readGeneOrder( geneOrderStr ):

	inFile = open( geneOrderStr, 'r' )
	geneAr = []
	
	hasVals = None
	
	for line in inFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split()
		
		if hasVals == None:
			hasVals = len(lineAr) > 1
		
		gName = lineAr[0]
		if hasVals:
			val = float(lineAr[1])
			geneAr += [(val, gName)]
		else:
			geneAr += [gName]
	# end for line
	
	inFile.close()
	if hasVals:
		geneAr.sort(reverse=True)
		return [ g[1] for g in geneAr ]
	else:
		return geneAr
		
def processFile( bamFileStr, sampleName, geneDict, geneOrder, distance, numBins, percentile, outId, isCompress, isScaled, isPrint, info ):
	isValid = checkBAM( bamFileStr )
	
	if not isValid:
		print('ERROR: Problem with {:s}..skipping'.format( os.path.basename(bamFileStr)) )
		return
	
	if isPrint:
		print('Begin processing',sampleName)
		
	outFileStr = outId + '_' if outId != None else ''
	outFileStr += sampleName + '_abs-tss_n' + str(numBins) + '_d' + str(distance) + '.tsv'
	outFileStr += '.gz' if isCompress else ''
	
	# otherwise open bam file
	bamFile = pysam.AlignmentFile( bamFileStr, 'rb' )
	scaleVal = float(bamFile.mapped / 1000000) if isScaled else 1.0
	
	outMatrix = []
	numList = []
	geneList = []
	
	# loop through genes
	for gene in geneOrder:
		geneInfo = geneDict.get( gene )
		if geneInfo == None:
			continue
		chrm, start, end, strand = geneInfo
		
		# only keep long enough genes
		if (end - start) < distance or start - distance < 0:
			continue
			
		tss = countRegion( bamFile, chrm, start, distance, numBins, scaleVal )
		tts = countRegion( bamFile, chrm, end, distance, numBins, scaleVal )
		
		# negative strand gene
		if strand == '-':
			tss.reverse()
			tts.reverse()
			outAr = tts + tss
		else:
			outAr = tss + tts
		
		outMatrix.append(outAr)
		numList += outAr
		geneList += [gene]
	# end for gene
	
	bamFile.close()
	
	# get percentile threshold
	threshold = calculatePercentile( numList, percentile )
	
	# write output
	if isPrint:
		print( 'Writing output for', sampleName, 'to', os.path.basename(outFileStr) )
		
	writeOutput( outFileStr, outMatrix, geneList, threshold, info, isCompress )
	
	return True

def checkBAM( bamFileStr ):
	bamIndex = bamFileStr + '.bai'
	if os.path.isfile(bamFileStr) == False:
		print('ERROR: BAM file does not exist')
		return False
	elif os.path.isfile(bamIndex) == False:
		print('WARNING: BAM index file does not exist...creating')
		pysam.index( bamFileStr )
	return True	
	
def countRegion( bamFile, chrm, midPoint, distance, numBins, scaleVal ):
		
	totalBins = numBins*2
	outAr = [0] * (totalBins)
	
	rStart = midPoint - distance
	#rEnd = midPoint + distance
	binSize = distance / numBins
		
	for i in range(totalBins):
		pStart = i*binSize + rStart
		pEnd = pStart + binSize
		reads = bamFile.fetch( chrm, pStart, pEnd )
		nReads = reduce(lambda x, y: x+1, reads, 0)
		outAr[i] = float(nReads) / scaleVal
	# end for i
	return outAr

def calculatePercentile( numList, percentile ):
	
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

def writeOutput( outFileStr, outMatrix, geneList, threshold, info, isCompress ):
	if isCompress:
		outFile = gzip.open( outFileStr, 'wt' )
	else:
		outFile = open( outFileStr, 'w' )
	
	headerAr = [ 'gene', 'gene.num', 'bin', 'value' ]
	outFile.write( info + '\n' + '\t'.join(headerAr) + '\n')
	
	nGenes = len(geneList)
	
	for i in range(nGenes):
		geneNum = 'G{:08}'.format( nGenes - i )
		tmpAr = outMatrix[i]
		
		for j in range(len(tmpAr)):
			outStr = '{:s}\t{:s}\t{:d}\t{:f}\n'.format( geneList[i], geneNum, j, (threshold if tmpAr[j] > threshold else tmpAr[j]) )
			outFile.write( outStr )
		# end for j
	# end for i
	outFile.close()

def parseInputs( argv ):
	isPrint = True
	isCompress = False
	isScaled = True
	numBins = NUMBINS
	distance = DISTANCE
	outId = None
	labels = None
	numProc = NUMPROC
	percentile = PERCENTILE
	geneOrderStr = None
	startInd = 0
	
	for i in range(len(argv)):
		if argv[i] == '-h':
			printHelp()
			exit()
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i] == '-z':
			isCompress = True
			startInd += 1
		elif argv[i] == '-n':
			isScaled = False
			startInd += 1
		elif argv[i].startswith( '-b=' ):
			try:
				numBins = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of bins must be an integer...using default', NUMBINS )
		elif argv[i].startswith( '-d=' ):
			try:
				distance = argv[i][3:]
				if distance.endswith( 'k' ) or distance.endswith( 'K' ):
					distance = int( distance[:-1] ) * 1000
				elif distance.endswith( 'm' ) or distance.endswith( 'M' ):
					distance = int( distance[:-1] ) * 1000000
				else:
					distance = int( distance )
				startInd += 1
			except ValueError:
				print( 'WARNING: distance must be an integer..using default', DISTANCE )
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be an integer...using default', NUMPROC )
		elif argv[i].startswith( '-t=' ):
			try:
				percentile = float( argv[i][3:] )
				if percentile > 1:
					percentile /= 100
				startInd += 1
			except ValueError:
				print( 'ERROR: percentile must be numeric' )
				exit()
		elif argv[i].startswith( '-r=' ):
			geneOrderStr = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outId = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-l='):
			labels = argv[3:].split(',')
			startInd += 1
		elif argv[i].startswith('-'):
			print( 'ERROR: unknown parameter {:s}. Use -h to see available parameters'.format( argv[i] ) )
			exit()
	# end for i
	
	gffFileStr = argv[startInd]
	bamFileAr = []
	for j in range( startInd + 1, len(argv) ):
		bamFileAr += [ argv[j] ]
	
	# check that numBins divisble
	if distance % numBins != 0:
		nnsize = int( math.ceil( float(distance) / numBins ) )
		distance = nnsize * numBins
		print( 'WARNING: num bins does not evenly divide into distance...adjusting distance to', distance )
	
	# check correct number of labels
	if labels != None and len(bamFileAr) != len(labels):
		print('WARNING: number of labels does not match number of BAM files. Ignoring labels parameter')
		labels = None
	
	# check if gene order file exists
	if os.path.isfile(geneOrderStr) == False:
		print( 'WARNING: gene order file {:s} does not exist...using default gene ordering'.format(os.path.basename(geneOrderStr)) )
		geneOrderStr = None
	
	processInputs( gffFileStr, bamFileAr, numProc, numBins, distance, outId, percentile, geneOrderStr, labels, isPrint, isCompress, isScaled )

def printHelp():
	print('Usage:\tpython heatmap_abs_bam_pe.py [-h] [-q] [-z] [-n] [-b=num_bins]')
	print('\t[-d=distance] [-o=out_id] [-t=percentile] [-r=gene_order] [-l=labels]')
	print('\t<gff_file> <bam_file> [bam_file]*')
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
