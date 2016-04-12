import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: python3 sequences_from_bed.py [-m] [-s=upstream[,downstream]] <fasta_file> <bed_file> [bed_file]*
# get references sequences corresponding to the regions in BED file
# outputs to a fasta file by default

STREAMSIZE=0
NUMPROC=1

def processInputs( fastaFileStr, bedFileAr, upstreamSize, downstreamSize, useMidpoint, numProc ):
	
	# read fasta
	print( 'Reading FASTA' )
	fastaFile = FileFASTA( fastaFileStr )
	fastaDict = fastaFile.readFasta()
	
	# go through bed files
	'''print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processBedFile, args=(f, fastaDict, upstreamSize, downstreamSize, useMidpoint ) ) for f in bedFileAr ]
	outs = [ p.get() for p in results ]'''
	for f in bedFileAr:
		 processBedFile( f, fastaDict, upstreamSize, downstreamSize, useMidpoint )
	print( 'Done' )

def processBedFile( bedFileStr, fastaDict, upstreamSize, downstreamSize, useMidpoint ):
	print( 'Reading', bedFileStr )
	outFileStr = getOutfileName( bedFileStr )
	bedFile = open( bedFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in bedFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		name = lineAr[3]
		chrm = lineAr[0]
		start = int( lineAr[1] ) + 1
		end = int( lineAr[2] ) + 1
		if fastaDict.get( chrm ) == None:
			print( 'WARNING: {:s} not included in fasta'.format( chrm ) )
			continue
		fastaSeq = fastaDict[chrm]
		# using midpoint
		if useMidpoint:
			mpoint = int( math.floor( float(end - start) / 2.0) ) + start
			sStart = max( 1, (mpoint - upstreamSize) )
			sEnd = min( len(fastaSeq), mpoint + downstreamSize )
			# check for not same length
			if sEnd - sStart != downstreamSize + upstreamSize:
				print( 'WARNING: {:s} not included because not same length'.format( name ) )
		# use whole region
		else:
			sStart = max(1, (start - upstreamSize))
			sEnd = min( len(fastaSeq), end + downstreamSize )
		#print( sStart, sEnd )
		outSeqAr = fastaSeq[ sStart:sEnd ]
		outStr = '>{:s} {:s} {:d} {:d}\n{:s}\n'.format( name, chrm, sStart, sEnd, ''.join(outSeqAr) )
		outFile.write( outStr )
	# end for line
	outFile.close()
	bedFile.close()
	print( 'Finished writing', outFileStr )
	return 1

def getOutfileName( bedFileStr ):
	
	bName = os.path.basename( bedFileStr )
	rInd = bName.rfind( '.' )
	if rInd != -1:
		bName = bName[:rInd]
	return bName + '_seq.fa'
	
def parseInputs( argv ):
	upstreamSize = STREAMSIZE
	downstreamSize = STREAMSIZE
	useMidpoint = False
	numProc = NUMPROC
	startInd = 0
	
	for i in range(min(3,len(argv))):
		if argv[i].startswith( '-s=' ):
			try:
				streamAr = argv[i][3:].split(',')
				upstreamSize = int( streamAr[0] )
				if len( streamAr ) == 2:
					downstreamSize = int( streamAr[1] )
				else:
					downstreamSize = int( streamAr[0] )
				startInd += 1
			except ValueError:
				print( 'ERROR: stream size must be integer' )
				exit()
		elif argv[i] == '-m':
			useMidpoint = True
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i] in ['-h', '--help', '-help' ]:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			eixt()
	# end for
	fastaFileStr = argv[startInd]
	bedFileAr = []
	
	for j in range(startInd+1, len(argv) ):
		bedFileAr += [ argv[j] ]
	
	processInputs( fastaFileStr, bedFileAr, upstreamSize, downstreamSize, useMidpoint, numProc )
	

def printHelp( ):
	print( 'Usage: python3 sequences_from_bed.py [-m] [-s=upstream[,downstream]] <fasta_file> <bed_file>' )
	

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
