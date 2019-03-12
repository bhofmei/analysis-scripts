import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python decode_single_pileup.py [-o=out_id] <pileup_file>

def processInputs( pileFileStr, outID ):
	print( 'Reading', os.path.basename( pileFileStr ) )
	outStr = parsePileup( pileFileStr )
	info = '#from_script: decode_single_pileup.py; pileup_file: ' + os.path.basename( pileFileStr )
	outFileStr = '{:s}_decoded_pileup.tsv'.format( outID )
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outStr, info )
	print( 'Done' )

def parsePileup( inFileStr ):
	inFile = open( inFileStr, 'r' )
	outStr = ''
	
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) pos (2) ref (3) read counts (4) nucs (5) qual
		refBase = lineAr[2]
		nucs = list(lineAr[4])
		nucDecode = [ decode( refBase, x ) for x in nucs ]
		nucCounts = countBases( nucDecode )
		outStr += '\t'.join( lineAr[0:4] ) + '\t' + '\t'.join ( [ str(x) for x in nucCounts[:-1] ] ) + '\n'
	# end for line
	inFile.close()
	return outStr
	
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

def writeOutput( outFileStr, outStr, info ):
	outFile = open( outFileStr, 'w' )
	headerAr = ['chrm','pos','refBase','readCounts','A', 'C', 'G', 'T', 'a', 'c', 'g', 't' ]
	outFile.write( info + '\n' + '\t'.join( headerAr ) + '\n' )
	outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	outID = 'out'
	startInd = 0
	
	for i in range(min(1,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	pileupFileStr = argv[startInd]
	processInputs( pileupFileStr, outID )


if __name__ == "__main__":
	if len(sys.argv) < 1 :
		print ("Usage: python decode_pileup.py [-o=out_id] <pileup_file>")
	else:
		parseInputs( sys.argv[1:] )
