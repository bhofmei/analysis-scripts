import sys, math, glob, multiprocessing, subprocess, os

# Usage: table_allc_pe.py [-cg] [-p=num_proc] [-c=chrm_list] <out_prefix> <allc_path> <sample1> <sample2> [sampleN]*
# creates table of samples and methylation per positions

NUMPROC=1
CHRMLIST=["Chr1","Chr2","Chr3","Chr4","Chr5"]

def processInputs( outPre, allcPath, sampleNamesAr, numProc, chrmList, isCG ):
	pool = multiprocessing.Pool( processes=numProc )
	# process chromosomes
	print( "Using", numProc, "processors" )
	results = [ pool.apply_async(processChrm, args=(chrm, outPre, allcPath, sampleNamesAr, isCG) ) for chrm in chrmList ]
	outs = [ p.get() for p in results ]
	print( "Done." )
	
def processChrm( chrm, outPre, allcPath, sampleNamesAr, isCG ):
	print( "Processing chrm", chrm )
	outFileStr = outPre + "_" + chrm 
	if isCG:
		outFileStr += "_cg"
	outFileStr += ".tsv"
	allcDict = {}
	
	# loop through samples
	for name in sampleNamesAr:
		allcFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, name, chrm ))
		print( "Reading", name, "-", chrm )
		allcDict = readAllc( allcFileStr, allcDict, isCG )
	# write output
	print( "Writing", outFileStr )
	writeOutput( outFileStr, allcDict, sampleNamesAr )

def readAllc( allcFileStr, allcDict, isCG ):
	
	allcFile = open( allcFileStr, 'r' )
	for line in allcFile:
		if line.startswith('c'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if (isCG and lineAr[3].startswith('CG')) or (isCG == False) :
			pos = int(lineAr[1])
			wMeth = float(lineAr[4]) / float(lineAr[5])
			if allcDict.get(pos)==None:
				allcDict[pos] = []
			allcDict[pos] += [ wMeth ]
	# end for line
	allcFile.close()
	return allcDict
	
def writeOutput( outFileStr, allcDict, sampleNamesAr ):
	
	n = len( sampleNamesAr )
	outFile = open( outFileStr, 'w' )
	header = "pos\t" + "\t".join(sampleNamesAr) + "\n"
	outFile.write( header )
	# loop through positions
	for pos in sorted( allcDict.keys() ):
		ar = allcDict[pos]
		# only include arrays that have enough values
		if len(ar) == n:
			arS = [ '{:.4f}'.format(x) for x in ar ]
			outStr = '{:d}\t{:s}\n'.format( pos, '\t'.join(arS) )
			outFile.write( outStr )
	# end for
	outFile.close()
		

def parseInputs( argv ):
	numProc = NUMPROC
	chrmList = CHRMLIST
	isCG = False
	startInd = 0
	
	for i in range(min(3,len(argv)-4)):
		if argv[i].startswith('-p='):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith('-c='):
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i] == '-cg':
			isCG = True
			startInd += 1
	# end for
	outPre = argv[startInd]
	allcPath = argv[startInd+1]
	sampleNamesAr = []
	for j in range(startInd+2, len(argv) ):
		sampleNamesAr += [ argv[j] ]
	processInputs( outPre, allcPath, sampleNamesAr, numProc, chrmList, isCG )

if __name__ == "__main__":
	if len(sys.argv) < 5 :
		print ("Usage: table_allc_pe.py [-cg] [-p=num_proc] <out_prefix> <allc_path> <sample1> <sample2> [sampleN]*")
	else:
		parseInputs( sys.argv[1:] )
