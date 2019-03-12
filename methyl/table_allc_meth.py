import sys, math, os

# Usage: table_allc_meth.py [-cg] [-p=num_proc] [-c=chrm_list] <out_prefix> <allc_path> <sample1> <sample2> [sampleN]*
# creates table of samples and methylation per positions

def processInputs( allcPath, sampleNamesFile, outID, minCov ):
	# read file
	sampleFilesAr, sampleNamesAr = readSampleFile( sampleNamesFile )
	print( 'Found {:d} samples'.format(len(sampleNamesAr)) )
	print( 'Min coverage: ', minCov )

	outFileStr = outID + ".tsv"
	allcDict = {}
	
	# loop through samples
	for file in sampleFilesAr:
		allcFileStr = os.path.normpath('{:s}/{:s}'.format( allcPath, file ))
		print( "Reading", file )
		allcDict = readAllc( allcFileStr, allcDict, minCov )
	
	# write output
	print( "Writing", outFileStr )
	writeOutput( outFileStr, allcDict, sampleNamesAr )
	print( 'Done' )
	
def readSampleFile( sampleFileStr ):
	
	outAr = []
	namesAr = []
	inFile = open( sampleFileStr, 'r' )
	for line in inFile:
		if line.startswith('#'):
			continue
		line = line.rstrip()
		outAr += [ line ]
		name = line.replace('allc_', '').replace('.tsv','')
		namesAr += [ name ]
	# end for
	inFile.close()
	return outAr, namesAr

def readAllc( allcFileStr, allcDict, minCov ):
	
	allcFile = open( allcFileStr, 'r' )
	for line in allcFile:
		lineAr = line.rstrip().split('\t')
		if len( lineAr ) < 7 or lineAr[6].isdigit() == False:
			continue
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		chrm = lineAr[0]
		# add chrm if not all ready added
		if allcDict.get(chrm) == None:
			allcDict[chrm] = {}
		# get position
		pos = int(lineAr[1])
		nreads = int( lineAr[5] )
		wMeth = int( lineAr[6] )
		if allcDict[chrm].get(pos)==None:
			allcDict[chrm][pos] = []
		if nreads >= minCov:
			allcDict[chrm][pos] += [ wMeth ]
	# end for line
	allcFile.close()
	return allcDict
	
def writeOutput( outFileStr, allcDict, sampleNamesAr ):
	
	n = len( sampleNamesAr )
	outFile = open( outFileStr, 'w' )
	header = "chrm\tpos\t" + "\t".join(sampleNamesAr) + "\n"
	outFile.write( header )
	# loop through chrms
	for chrm in sorted( allcDict.keys() ):
		chrmDict = allcDict[chrm]
		for pos in sorted( chrmDict.keys() ):
			ar = chrmDict[pos]
			# only include arrays that have enough values
			if len(ar) == n:
				arS = [ '{:d}'.format(x) for x in ar ]
				outStr = '{:s}\t{:d}\t{:s}\n'.format( chrm, pos, '\t'.join(arS) )
				outFile.write( outStr )
	# end for
	outFile.close()	

def parseInputs( argv ):
	startInd = 0
	minCov = 1
	outID = 'out'
	
	for i in range(min(2,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum coverage must be integer' )
				exit()
	# end for
	
	#outPre = argv[startInd]
	allcPath = argv[startInd]
	sampleNamesFile = argv[startInd+1]
	processInputs( allcPath, sampleNamesFile, outID, minCov )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: table_allc_meth.py [-o=out_id] [-v=min_cov] <allc_path> <file_of_samples>\nout_id\tprefix for output file [default out]\nmin_cov\tmin read cov to be include [default 1]\nallc_path\tpath to allc files\nfile_of_samples\ttext file with each allc file listed per line")
	else:
		parseInputs( sys.argv[1:] )
