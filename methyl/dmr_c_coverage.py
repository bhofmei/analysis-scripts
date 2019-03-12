import sys, math, glob, multiprocessing, subprocess, os

# Usage: dmr_c_coverage.py [-p=num_proc] [-o=out_prefix] <dmr_file> <fasta_file> <allc_path> <sample> [sampleN]*

NUMPROC=1

def processInputs( allcPath, sampleNamesAr, dmrFileStr, fastaFileStr, outPre, numProc ):
	print( 'reading dmr file' )
	dmrDict, dmrCount = readDMRFile( dmrFileStr )
	print( 'found', str(dmrCount), 'dmrs' )
	# number of c and g per region using fasta
	print( 'reading and processing fasta file' )
	fastaDict = readFasta( fastaFileStr )
	dmrArray, dmrStr = processDMRFasta( dmrDict, fastaDict )
	#print( dmrArray )
	del( fastaDict )
	
	# number of c per region using allc
	print( 'processing with', str(numProc), 'processors' )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(processSample, args=(name, allcPath, dmrDict)) for name in sampleNamesAr ]
	resultAr = [ p.get() for p in results ]
	#for i in range(len(resultAr[0])):
		#print( i, dmrStr[i], resultAr[0][i] )
	
	# write output
	outFileStr = outPre + '_c_coverage.tsv'
	print( 'writing output' )
	writeOutput( outFileStr, resultAr, sampleNamesAr, dmrArray, dmrStr )

def readDMRFile( dmrFileStr ):
	'''
		return dictionary of subset (DMR) regions
		{chr1:[(start,end),(start2,end2)],chr2:[(start,end)],...}
	'''
	dmrFile = open( dmrFileStr, 'r' )
	dmrDict = {}
	dmrCount = 0
	
	for line in dmrFile:
			lineAr = line.rstrip().split()
			chrm = lineAr[0]
			start = int( lineAr[1] )
			end = int( lineAr[2] )
			if dmrDict.get(chrm) == None:
				dmrDict[chrm] = []
			dmrDict[chrm] += [(start, end)]
			dmrCount += 1
	dmrFile.close()
	return dmrDict, dmrCount

def readFasta( fastaFileStr ):
	''' 
		we will read in the whole fasta because it shouldn't be too big
		read into dictionary where chrm is the key and and the value is the
		sequence
		NOTE: we are making the sequence 1-based indexed and automatically capitalizing
		all of the bases
		returns the created dictionary
	'''
	fastaFile = open( fastaFileStr, 'r' )
	fastaDict = {}
	chrm = None
	seq = [' ']
	for line in fastaFile:
		line = line.rstrip()
		# chromosome headers
		if line.startswith('>'):
			# need to write old sequence
			if seq != [' '] and chrm != None:
				fastaDict[chrm] = seq
				seq = [' ']
			lineAr = line.split(' ')
			chrm = lineAr[0].replace('>', '')
		
		# sequence
		else:
			seqL = list(line.upper())
			seq += seqL
			
	# handle last chrm read
	if chrm not in fastaDict:
		fastaDict[chrm] = seq
	fastaFile.close()
	return fastaDict

def processDMRFasta( dmrDict, fastaDict ):
	dmrArray = []
	dmrStr = []
	dictKeys = sorted( dmrDict.keys() )
	
	# loop through chrms
	for chrm in dictKeys:
		regionAr = dmrDict[chrm]
		fasta = fastaDict[chrm]
		for region in regionAr:
			seq = fasta[ region[0]: region[1]+1 ]
			cCG = seq.count('C') + seq.count( 'G' )
			dmrArray += [ cCG ]
			dmrStr += [ '{:s}:{:d}-{:d}'.format( chrm, region[0], region[1] ) ]
	return dmrArray, dmrStr
	
def processSample( name, allcPath, dmrDict ):
	
	outAr = []
	# loop through chromosome
	dictKeys = sorted( dmrDict )
	for chrm in dictKeys:
		allcFileStr = allCFileStr =os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, name, chrm ) )
		# read allc
		print( 'reading and processing allc for', name, chrm )
		allcDict = readAllc( allcFileStr )
		# loop through regions
		for region in dmrDict[chrm]:
			r, c = countC( allcDict, region[0], region[1] )
			outAr += [ (r,c) ]
	# end for chrm
	return outAr

def readAllc( allcFileStr ):
	allcDict = {}
	allcFile = open( allcFileStr, 'r' )
	header = True
	for line in allcFile:
		if header:
			header = False
			continue;
		lineAr = line.rstrip().split( '\t' )
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		pos = int( lineAr[1] )
		allcDict[pos] = int( lineAr[5] )
	# end for
	allcFile.close()
	return allcDict
	
def countC( allcDict, start, end ):
	
	count = 0
	reads = 0
	for i in range(start, end+1):
		t = allcDict.get(i)
		if t != None:
			reads += t
			count += 1
	return reads, count

def writeOutput( outFileStr, resultAr, sampleNamesAr, dmrArray, dmrStr ):
	# dmr region sample c_reads total_c c_cov
	header = 'DMR\tregion\tsample\tc.with.read\tc.reads\tc.total\tc.cov\n'
	outFile = open( outFileStr, 'w' )
	outFile.write( header )
	
	# loop through sample
	for i in range(len(sampleNamesAr)):
		# loop through region
		for j in range(len(dmrStr)):
			cov = float( resultAr[i][j][0] ) / dmrArray[j]
			outStr = '{:d}\t{:s}\t{:s}\t{:d}\t{:d}\t{:d}\t{:.4f}\n'.format( j, dmrStr[j], sampleNamesAr[i], resultAr[i][j][1], resultAr[i][j][0], dmrArray[j], cov )
			outFile.write( outStr )
		# end for j
	# end for i
	outFile.close()
	
def parseInputs( argv ):
	numProc = NUMPROC
	outPre = 'out'
	startInd = 0
	
	for i in range(min(2,len(argv)-4)):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith('-o='):
			outPre = argv[i][3:]
			startInd += 1
	# end for
	
	dmrFileStr = argv[startInd]
	fastaFileStr = argv[startInd + 1]
	allcPath = argv[startInd + 2]
	sampleNamesAr = []
	for j in range(startInd + 3, len(argv)):
		sampleNamesAr += [ argv[j] ]
	processInputs( allcPath, sampleNamesAr, dmrFileStr, fastaFileStr, outPre, numProc )

if __name__ == "__main__":
	if len(sys.argv) <5 :
		print ("Usage: dmr_c_coverage.py [-p=num_proc] [-o=out_prefix] <dmr_file> <fasta_file> <allc_path> <sample> [sampleN]*")
	else:
		parseInputs( sys.argv[1:] )
