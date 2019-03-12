import sys, multiprocessing, subprocess, os

# Usage: python decode_pileup_pe.py [-o=out_id] [-p=num_proc] <pileup_file> [pileup_file]*
## Assumes all pileup files have information for the same positions

NUMPROC=1

def processInputs( pileFileStrAr, outID, numProc, isPrint ):
	baseNamesAr = [ os.path.basename( x ) for x in pileFileStrAr ]
	info = '#from_script: decode_pileup_pe.py; pileup_files: ' + ','.join( baseNamesAr )
	
	if isPrint:
		print( 'Input files:', ', '.join( baseNamesAr ) )
	# get position information from first file
	posMat = parsePileupPos( pileFileStrAr[0] )
	nPos = len( posMat )
	
	# get sample names
	sampleNamesAr = getSampleNames( pileFileStrAr )
	
	if isPrint:
		print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(x, nPos) ) for x in pileFileStrAr ]
	outMat = [ p.get() for p in results ]
	# need to transpose
	outMat2 = transpose( outMat )
	
	if len(outMat2) != nPos and isPrint:
		print( 'WARNING: npos and outmat dimension issue' )
		print( nPos )
		print( len(outMat2), len(outMat2[0]) )
	
	outFileStr = '{:s}_decoded_pileup.tsv'.format( outID )
	if isPrint:
		print( 'Writing output to', outFileStr)
	writeOutput( outFileStr, posMat, outMat2, sampleNamesAr, info )
	if isPrint:
		print( 'Done' )

def parsePileupPos( inFileStr ):
	outMat = []
	inFile = open( inFileStr, 'r' )
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) pos (2) ref (3) read counts (4) nucs (5) qual
		outMat.append( lineAr[0:3] )
	# end for line
	inFile.close()
	return outMat

def getSampleNames( pileFileStrAr ):
	sampleNamesAr = []
	
	for file in pileFileStrAr:
		bName = os.path.basename( file )
		rInd = bName.rfind( '.' )
		if rInd != -1:
			bName = bName[:rInd]
		bName = bName.replace('_pileup','').replace('-pileup','')
		sampleNamesAr += [ bName ]
	# end for
	return sampleNamesAr

def processFile( pileFileStr, nPos ):
	print( 'Reading', os.path.basename( pileFileStr ) )
	outAr = parsePileup( pileFileStr )
	if len( outAr ) != nPos:
		print( 'WARNING: file {:s} has a different number of positions listed...using first {:d}'.format( os.path.basename(pileFileStr), nPos ) )
		if len(outAr) > nPos:
			outAr = outAr[:nPos]
		else:
			outAr += [ 'NA' ] * ( nPos - len(outAr) )
	return outAr

def parsePileup( inFileStr ):
	inFile = open( inFileStr, 'r' )
	outAr = []
	
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) pos (2) ref (3) read counts (4) nucs (5) qual
		refBase = lineAr[2]
		nucs = list(lineAr[4])
		nucDecode = [ decode( refBase, x ) for x in nucs ]
		nucCounts = countBases( nucDecode )
		outAr += [ formatBasesOutput( nucCounts ) ]
	# end for line
	inFile.close()
	return outAr
	
def decode( refBase, n ):
	baseAr = ['A', 'C', 'G', 'T' ]
	revAr = [ 't', 'g', 'c', 'a' ]
	anyBase = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't' ]
	
	try:
		ind = baseAr.index( refBase )
		if n == '.':
			return baseAr[ind]
		elif n == ',':
			return revAr[ind]
		elif n in anyBase:
			return n
		else:
			return ''
	except ValueError:
		return 'X'

def countBases( inAr ):
	outAr = [0] * 9
	anyBase = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'X' ]
	for n in inAr:
		if n == '':
			continue
		ind = anyBase.index( n )
		outAr[ind] += 1
	# end for n
	return outAr

def formatBasesOutput( inAr ):
	anyBase = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't' ]
	outStr = ''
	for i in range( 8 ):
		c = inAr[i]
		if c != 0:
			outStr += '{:s}({:d}) '.format( anyBase[i], inAr[i] )
	# end for i
	# remove last space
	return outStr[:-1]

def transpose( inMat ):
	nrow = len( inMat )
	ncol = len( inMat[0] )
	outMat = [ [None] * nrow for x in range(ncol) ]
	for i in range(nrow):
		for j in range(ncol):
			outMat[j][i] = inMat[i][j]
		# end for j
	# end for i
	return outMat

def writeOutput( outFileStr, posMat, outMat2, sampleNamesAr, info ):
	outFile = open( outFileStr, 'w' )
	headerAr = ['chrm','pos','refBase'] + sampleNamesAr
	outFile.write( info + '\n' + '\t'.join( headerAr ) + '\n' )
	for i in range( len(posMat) ):
		outStr = '\t'.join( posMat[i] )
		outStr += '\t' + '\t'.join( outMat2[i] )
		outStr += '\n'
		outFile.write( outStr )
	# end for i
		
	outFile.close()

def parseInputs( argv ):
	outID = 'out'
	numProc = NUMPROC
	isPrint = True
	startInd = 0
	
	for i in range(min(3,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				
			except ValueError:
				print( 'WARNING: number of processors must be integer...using', NUMPROC )
				numProc = NUMPROC
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	pileupFileStrAr = []
	for j in range( startInd, len(argv) ):
		pileupFileStrAr += [argv[j]]
	
	processInputs( pileupFileStrAr, outID, numProc, isPrint )

def printHelp():

	print( 'Usage:\tpython decode_pileup_pe.py [-h] [-q] [-o=out_id] [-p=num_proc]\n\t<pileup_file> [pileup_file]*' )
	print()
	print( 'Required:' )
	print( 'pileup_file\tpileup file for a sample; output from samtools pileup' )
	print()
	print( 'Optional:' )
	print( '-h\t\tprint help and exit' )
	print( '-q\t\tquiet; do not print progress' )
	print('-o=out_id\tidentifier for output file [default "out"]' )
	print( '-p=num_proc\tnumber of processors [default {:d}]'.format( NUMPROC) )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
