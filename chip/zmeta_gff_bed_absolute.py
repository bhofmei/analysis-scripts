import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 meta_gff_bed_absolute.py <gff_file> <distance> <bed_file> [bed_file]*

BINSIZE=100

def readGFF( gffFileStr ):
	gffFile = open (gffFileStr, 'r' )
	gffAr = []
	
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
		
		# change this to select a different gff feature 
		# or remove and un-indent line after
		if lineAr[2] == "gene":
			gffAr+= [(scaffold, start, end, strand)]
	
	gffFile.close()
	return gffAr

def readBed(bedFileStr ):
	''' takes in the bed file, scaffold we are considering
		returns a dictionary where the position is the key and the value is the number
		of hits for that position
	'''
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	readCount = 0
	
	# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
	for line in bedFile:
		line = line.replace('\n','')
		lineAr =line.split('\t')
		try:
			curScaf = lineAr[0]
			pos = int(lineAr[1])
		except ValueError:
			pass
		else:
			readCount += 1
			# no dictionary for scaffold
			if bedDict.get(curScaf) == None:
				bedDict[curScaf] = {pos:1}
			# dictionary for scaffold but position not included
			elif bedDict.get(curScaf).get(pos) == None:
				bedDict[curScaf][pos] = 1
			# dictionary for scaffold and position included
			else:
				bedDict[curScaf][pos] += 1
			#print bedDict[curScaf]
	
	bedFile.close()
	return bedDict, readCount

def processBed( bedDict, gffAr, distance, readCount ):
	#print( bedDict['Chr13'] )
	outArF = [ [0,0] for x in range( int( distance // BINSIZE ) * 2) ]
	outArT = [ [0,0] for x in range( int( distance // BINSIZE ) * 2) ]
	
	for gene in gffAr:
		scaffold = gene[0]
		bedDictS = bedDict.get( scaffold )
		if bedDictS == None:
			continue
		start = gene[1]
		end = gene[2]
		strandNum = ( 0 if gene[3] == '+' else 1 )
		curArF, curArT = bedCountRegion( bedDictS, start, end, distance, strandNum )
		outArF = addToArray( outArF, curArF )
		outArT = addToArray( outArT, curArT )
	outAr = correctArray( outArF, readCount ) + [0] * 2 + correctArray( outArT, readCount )	
	return outAr

def bedCountRegion( bedDict, start, end, distance, strandNum ):
	numBins = int(distance//BINSIZE)
	outArF = [-1] * (numBins*2)
	outArT = [-1] * (numBins*2)
	geneLength = end - start + 1
	
	# long enough genes
	if geneLength >= 2*distance:
		for bin in range( len( outArF ) ):
			for i in range( BINSIZE ):
				# 5` for + genes
				if outArF[bin] == -1:
					outArF[bin] = 0
				dictEntry = bedDict.get((start - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArF[bin] += dictEntry
				# 3` for + genes
				if outArT[bin] == -1:
					outArT[bin] = 0
				dictEntry = bedDict.get((end - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArT[bin] += dictEntry
	else:
		indStop = int((geneLength / 2 ) // BINSIZE)
		ind1 = 0
		ind2 = numBins  + indStop
		ind3 = numBins - indStop
		ind4 = len( outArT )
		
			
		# 5` for + genes (3` for - genes)
		for bin in range( ind1, ind2 ):
			for i in range( BINSIZE ):
				if outArF[bin] == -1:
					outArF[bin] = 0
				dictEntry = bedDict.get((start - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArF[bin] += dictEntry
		# 3` for + genes (5` for - genes)
		for bin in range( ind3, ind4 ):
			for i in range( BINSIZE ):
				if outArT[bin] == -1:
					outArT[bin] = 0
				dictEntry = bedDict.get((end - distance) + (bin*numBins) + i)
				if dictEntry != None:
					outArT[bin] += dictEntry
					
	# need to reverse for negative strand and return opposite
	if strandNum == 1:
		outArF.reverse()
		outArT.reverse()
		return outArT, outArF
	else:
		return outArF, outArT

def writeOutput( outFileStr, outMatrix, sampleNames ):
	
	outFile = open( outFileStr, 'w' )
	
	# headers would be: bin, count, type
	for j in range( len(outMatrix) ):
		for i in range( len(outMatrix[j]) ):
			outStr = '{:d},{:f},{:s}\n'.format(i+1, outMatrix[j][i], sampleNames[j] )
			outFile.write( outStr )
	outFile.close()

def addToArray( oldArray, newArray ):
	'''
		add to meta array
		oldArray is an array of 2 element arrays [bedCount, numGenes]
		newArray is array of bedCount
	'''
	if len( oldArray ) != len( newArray ):
		return None
	else:
		#print ( newArray )
		for i in range( len( newArray ) ):
			if newArray[i] != -1:
				oldArray[i][0] += newArray[i]
				oldArray[i][1] += 1
	return oldArray

def correctArray( inAr, readCount ):
	outAr = []
	
	for i in range( len(inAr) ):
		outAr += [ (float(inAr[i][0])*1000) / (inAr[i][1]*readCount) ]
	return outAr

def getSampleName( fileStr ):
	lind = fileStr.rfind('/')
	rind = fileStr.rfind('.')
	name = fileStr[lind+1:rind]
	return name
	
def processInputs( gffFileStr, distance, bedFileStrAr ):
	
	gffAr = readGFF( gffFileStr )
	print( 'read gff file' )
	
	sampleNames = [ getSampleName( f ) for f in bedFileStrAr ]
	outMatrix = []
	
	for bedFileStr in bedFileStrAr:
		bedDict, readCount = readBed( bedFileStr )
		print( 'read', bedFileStr )
		metaAr = processBed( bedDict, gffAr, distance, readCount )
		outMatrix.append( metaAr )
	
	writeOutput( 'meta_absolute_'+str(distance)+'.csv', outMatrix, sampleNames )
	print( 'output written' )
	print( 'Done.' )

def parseInputs( argv ):
	gffFileStr = argv[0]
	try:
		distance = int( argv[1] )
		if distance % BINSIZE != 0:
			print( 'ERROR: distance must be divisible by', BINSIZE )
			exit()
	except ValueError:
		print( 'ERROR: distance must be an integer' )
		exit()
	bedFileStrAr = []
	for i in range( 2, len( argv ) ):
		bedFileStrAr += [ argv[i] ]
	processInputs( gffFileStr, distance, bedFileStrAr )
	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 meta_gff_bed_absolute.py <gff_file> <distance> <bed_file> [bed_file]*")
	else:
		parseInputs( sys.argv[1:] )
