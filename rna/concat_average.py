import sys, math, os, glob

# Usage: concat_average.py <outFile prefix> <fpkmFile1> [fpkmFileN]*

# keep count of number of sample on
# keep array of sample names
# read in first file - create general dictionary
# for each subsequent file, create individual dictionary
	# iterate through general dictionary in indiv dictionary
		# if there's a match, add fpkm for gene
			# remove dictionary entry
		# if there's not a match, add none or empty string
	# iterate through remaining genes in the individual dictionary
		# add to general dictionary - start array with enough empties
# finally, average for each gene and add to dictionary for each gene
# print to output

def generateFPKMDict( fpkmFileAr ):

	# dictionary for first file
	fpkmDict = readFPKMFile( fpkmFileAr[0] )
	sampleNum = 1
	#print ( printDict( fpkmDict ) )
	# add to dictionary for rest of files
	while sampleNum < len( fpkmFileAr ):
		tempDict = readFPKMFile( fpkmFileAr[ sampleNum ] )
		#print( 'SAMPLE NUM ', sampleNum )
		#print ( printDict( tempDict ) )
		fpkmDict = addToFPKMDict( fpkmDict, tempDict, sampleNum )
		#print( 'JOINED DICTIONARY' )
		#print ( printDict( fpkmDict ) )
		sampleNum += 1
		
	return fpkmDict

def readFPKMFile( fpkmFileStr ):
	indivDict = {}
	
	fpkmFile = open( fpkmFileStr, 'r' )
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split('\t')
	
		# (0) tracking_id (1) class_code (2) nearest_ref (3) gene_id (4) gene_short_id
		# (5) tss_id (6) locus (7) length (8) coverage (9) FPKM (10) FPKM_conf_low
		# (11) FPKM_conf_high (12) FPKM_status
		
		if lineAr[12] == 'OK' and lineAr[4] != '-':
			geneAr = lineAr[4].split(',')
			for x in geneAr:
				indivDict[ x ] = [float( lineAr[9] )]
	
	fpkmFile.close()
	return indivDict
	
def addToFPKMDict( fpkmDict, indivDict, sampleNum ):

	# iterate through fpkmDict
	for key in fpkmDict.keys():
		t = indivDict.get(key)
		
		# no match
		if t == None:
			fpkmDict[key] += [None]
		else:
			fpkmDict[key] += t
			del indivDict[key]
	
	# iterate through the rest of indivDict
	for key in indivDict.keys():
		fpkmAr = [None] * sampleNum
		fpkmAr += indivDict[key]
		fpkmDict[key] = fpkmAr

	return fpkmDict

def calculateAverage( fpkmDict ):
	
	for key in fpkmDict.keys():
		sum = 0
		numN = 0
		fpkmAr = fpkmDict[key]
		
		for i in fpkmAr:
			if i != None:
				sum += i
				numN += 1
		ave = float( sum / float(numN) )
		fpkmDict[key] += [ ave ]
		
	return fpkmDict

def calculateStandardDeviation( fpkmDict ):
	
	for key in fpkmDict.keys():
		sumSq = 0
		numN = 0
		fpkmAr = fpkmDict[key]
		# array is [fpkm1, ..., fpkmN, average]
		av = fpkmAr[-1]
		
		for i in fpkmAr[0:-1]:
			if i != None:
				sumSq += pow( (i - av), 2 )
				numN += 1
		sd = math.sqrt( float( sumSq / float( numN-1 ) ) )
		fpkmDict[key] += [sd]
	return fpkmDict
		

def printDict( dict ):
	keys = list(sorted(dict.keys()))
	size = len( keys )
	
	s = '-' * 50 + '\n'
	s += "size: {:d}\n".format(size)
	for i in range( 20 ):
		s += "{:s}:\t{:s}\n".format( keys[i] ,str( dict[keys[i]] ) )
	s += '-' * 50 + '\n'
	return s
	
	

def writeOutput ( outFileStr, fpkmDict, sampleAr ):
	outFile = open( outFileStr, 'w' )
	
	# write header
	s = 'gene,'
	for i in sampleAr:
		s += i + ','
	s += 'average' + '\n'
	outFile.write(s)
	
	# loop through dictionary
	for key in sorted( fpkmDict.keys() ):
		s = key
		fpkmAr = fpkmDict[key]
		
		for i in fpkmAr:
			if i == None:
				s+= ',NA'
			else:
				s += "," + "{:.4f}".format(i)
		s += '\n'
		
		outFile.write(s)
	
	outFile.close()
	
def parseFPKMFiles( fpkmFileAr ):
	sampleAr = []
	
	# check for directory 
	if os.path.isdir( fpkmFileAr[0] ):
			fpkmDir = fpkmFileAr[0]
			fpkmFileAr = glob.glob( fpkmDir + '/*.fpkm_tracking' )
	
	for s in fpkmFileAr:
		rind = s.rindex( '.' )
		s2 = s[0:rind]
		rind2 = s2.rindex( 'S' )
		sampleAr += [ s2[rind2:] ]
		
	return fpkmFileAr, sampleAr

def run( outFilePre, fpkmFileAr ):

	outFileStr = outFilePre + '.csv'
	fpkmFileAr, sampleAr = parseFPKMFiles( fpkmFileAr )
	print( 'Samples included: {:s}'.format( ' '.join( sampleAr ) ) )
	fpkmDict = generateFPKMDict( fpkmFileAr )
	fpkmDict = calculateAverage( fpkmDict )
	writeOutput ( outFileStr, fpkmDict, sampleAr )
	

if __name__ == "__main__":
	if len( sys.argv ) < 3:
		print( "python3.4 concat_average.py <outFile prefix> <fpkmFile1> [fpkmFileN]*\nor\npython3.4 concat_average.py <outFile_prefix> <fpkm_folder>" )
	else:
		outFilePre = sys.argv[1]
		fpkmStrAr = sys.argv[2:]
		run( outFilePre, fpkmStrAr )
		
	
