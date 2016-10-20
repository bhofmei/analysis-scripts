import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python pileup_genotype_pe.py [-o=out_id] [-p=num_proc] [-m=mother_label] [-f=father_label] [-v=min_cov] <decoded_pileup_file> 
## Assumes all pileup files have information for the same positions

NUMPROC=1
NUCAR = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't' ]

def processInputs( pileFileStr, outID, numProc, parentLabelAr, minCov ):
	
	info = '#from_script: pileup_genotype_pe.py; mother_label: {:s}; father_label: {:s}; min_cov: {:d}; input_files: {:s}'.format( parentLabelAr[0], parentLabelAr[1], minCov, os.path.basename(pileFileStr) )
	
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	print( 'Min coverage:', minCov )
	print( 'Input file:', os.path.basename(pileFileStr) )
	
	# parse input file
	posList, sampleNamesAr, valueMat = parseInputFile( pileFileStr )
	
	# check parents
	if parentLabelAr[0] not in sampleNamesAr or parentLabelAr[1] not in sampleNamesAr:
		print( 'ERROR: mother and/or father samples not found...use label parameters' )
		exit()
	parentInd = [ sampleNamesAr.index(x) for x in parentLabelAr ]
	
	# analyze valueMat with processors -> each row
	
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processRow, args=(valueMat[i], parentInd, minCov) ) for i in range(len(valueMat)) ]
	outMat = [ p.get() for p in results ]
	
	if outID == None:
		bName = os.path.basename( pileFileStr )
		rInd = bName.rfind( '.' )
		if rInd != -1:
			outID = bName[:rInd]
		else:
			outID = bName
	outFileStr = '{:s}_genotyped.tsv'.format( outID )
	print( 'Writing output to', outFileStr)
	writeOutput( outFileStr, posList, outMat, sampleNamesAr, info )
	print( 'Done' )
	
def parseInputFile( pileFileStr ):
	posList = []
	sampleNamesAr = []
	outMat = []
	sampleNamesColInt = 3
	isHeader = True
	inFile = open( pileFileStr, 'r' )
	
	for line in inFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split('\t')
		if isHeader:
			sampleNamesAr = lineAr[ sampleNamesColInt: ]
			isHeader = False
		else:
			# get position
			posList.append( lineAr[0:2] )
			outMat.append( lineAr[ sampleNamesColInt: ] )
	# end for line
	inFile.close()
	return posList, sampleNamesAr, outMat

def processRow( inAr, parInd, minCov ):
	n = len(inAr)
	outAr = [None] * n
	
	# get parent information
	mVal = decode( inAr[parInd[0]], minCov )
	fVal = decode( inAr[parInd[1]], minCov )
	#print( mVal, '--', fVal )
	pDiff = bitWise( mVal, fVal, 'xor' )
	# if no differences, return
	if sum(pDiff) == 0:
		return None
	# otherwise iterate through all samples and assign genotype
	mAr = bitWise( mVal, pDiff, 'and' )
	fAr = bitWise( fVal, pDiff, 'and' )
	
	for i in range(n):
		decoded = assignGenotype( inAr[i], minCov, mAr, fAr )
		if decoded == None:
			return None
		outAr[i] = ( inAr[i], decoded )
	# end for i
	return outAr
	
def bitWise( inAr1, inAr2, command ):
	outAr=[]
	if command == 'xor':
		outAr= [ inAr1[i] ^ inAr2[i] for i in range(len(inAr2)) ]
	elif command == 'and':
		outAr= [ inAr1[i] & inAr2[i] for i in range(len(inAr2)) ]
	elif command == 'or':
		outAr= [ inAr1[i] | inAr2[i] for i in range(len(inAr2)) ]
	return outAr
	
def decode( inStr, minCov ):
	outAr = [0] * 8
	anyBase = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']
	# split string by space
	strAr = inStr.split( ' ' )
	for x in strAr:
		# is of format nuc(count)
		fInd = anyBase.index( x[0] )
		c = int( x[x.find("(")+1:x.find(")")] )
		#print( fInd, c )
		outAr[fInd] = ( 1 if int(c) > minCov else 0 )
	# end for
	return outAr

def assignGenotype( sampleStr, minCov, mAr, fAr ):
	# decode new sample
	sampleAr = decode( sampleStr, minCov )
	#isM = sum( bitWise( sampleAr, mAr, 'and' ) ) == sum( mAr )
	#tm = 
	isM = sum( bitWise( sampleAr, mAr, 'and' )) > 0
	#isF = sum( bitWise( sampleAr, fAr, 'and' ) ) == sum( fAr )
	isF = sum(bitWise(sampleAr, fAr, 'and')) > 0
	
	if isM and isF:
		return 'MPV'
	elif isM:
		return 'mother'
	elif isF:
		return 'father'
	else:
		return None

def writeOutput( outFileStr, posList, outMat, sampleNamesAr, info  ):
	outFile = open( outFileStr, 'w' )
	headerAr = ['chrm','pos','sample', 'pileup', 'genotype']
	outFile.write( info + '\n' + '\t'.join( headerAr ) + '\n' )
	
	for i in range( len(posList) ):
		# get outMat values
		valAr = outMat[i]
		if valAr == None:
			continue	# skip position as it isn't useful
		# loop through samples
		for j in range(len(sampleNamesAr)):
			outStr = '\t'.join( posList[i] ) + '\t' + sampleNamesAr[j] + '\t' + '\t'.join( valAr[j] ) + '\n'
			outFile.write( outStr )
		# end for j
	# end for i
		
	outFile.close()

def parseInputs( argv ):
	outID = None
	numProc = NUMPROC
	parentLabelAr = ['mother', 'father']
	startInd = 0
	minCov = 0
	
	for i in range(min(5,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			parentLabelAr[0] = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			parentLabelAr[1] = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be integer...using 1' )
				numProc = NUMPROC
		elif argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: min coverage must be integer...using 0' )
				minCov = 0
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	pileupFileStr = argv[startInd]
	processInputs( pileupFileStr, outID, numProc, parentLabelAr, minCov )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python pileup_genotype_pe.py [-o=out_id] [-p=num_proc] [-m=mother_label] [-f=father_label] <decoded_pileup_file> ")
	else:
		parseInputs( sys.argv[1:] )
