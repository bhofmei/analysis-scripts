import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *
import bth_util

BINSIZE=500

# Usage: python3 fasta_nuc_density.py [-b=bin_size] [-o=outID] <nuc_search> <fasta_file>

def processInputs( fastaFileStr, nSearch, binSize, outID ):
	
	# check nSearch for only acgt
	nSearch = nSearch.upper()
	search = checkSearchString( nSearch )
	if  search == False:
		print( 'ERROR: search string invalid' )
		exit()
	print( 'Searching for {:s}'.format( ', '.join(search) ) )
	# read fasta
	print( 'Reading FASTA' )
	fastaFile = FileFASTA( fastaFileStr )
	fastaDict = fastaFile.getFastaDict()
	
	# loop through fasta chromosomes
	totalOutAr = []
	print( 'Searching chromosomes' )
	for chrm in sorted( fastaDict.keys() ):
		print( '-' + chrm )
		totalOutAr +=  searchFasta( fastaDict[chrm], search, binSize, chrm )
	
	# write output
	print( 'Writing output' )
	outFileStr = outID + '_' + nSearch + '_' + str(binSize) + '.tsv'
	writeOutput( outFileStr, totalOutAr )
	print( 'Done' )

def checkSearchString( input ):
	outAr = [input]
	for x in input:
		if x not in ['A', 'C', 'G', 'T', 'U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N' ]:
			return False
	# handle combinations
	if 'U' in input:
		outAr += transform( outAr,'U', ['T'] )
	if 'R' in input:
		outAr += transform( outAr, 'R', ['A','G'] )
	if 'Y' in input:
		outAr += transform( outAr,'Y', ['C','T'] )
	if 'S' in input:
		outAr += transform( outAr,'S', ['C','G'] )
	if 'W' in input:
		outAr += transform( outAr,'W', ['A','T'] )
	if 'K' in input:
		outAr += transform( outAr,'K', ['T','G'] )
	if 'M' in input:
		outAr += transform( outAr,'M', ['A','C'] )
	if 'B' in input:
		outAr += transform( outAr,'B', ['C','G','T'] )
	if 'D' in input:
		outAr += transform( outAr, 'D', ['A','G','T'] )
	if 'H' in input:
		outAr += transform( outAr, 'H', ['A','C','T'] )
	if 'V' in input:
		outAr += transform( outAr, 'V', ['A','G','C'] )
	if 'N' in input:
		outAr += transform( outAr, 'N', ['A','G','T','C'] )
	
	# clean
	return clean( outAr )


def transform(listAr, nucCode, options):
	outAr = []
	testAr = listAr
	while len(testAr) > 0:
		#print(testAr)
		if nucCode not in testAr[0]:
			if testAr[0] not in outAr:
				outAr += [testAr[0]]
			testAr = testAr[1:]
		else:
			testing = testAr[0]
			testAr = testAr[1:]
			for i in range(len(testing)):
				if testing[i] == nucCode:
					testAr += [ testing[:i] + x + testing[i+1:] for x in options ]
	return outAr

def clean( inAr ):
	outAr = []
	letters=['U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
	for x in inAr:
		if any(l in x for l in letters):
			continue
		else:
			outAr += [x]
	return outAr

def searchFasta( inAr, search, binSize, chrmName ):
	
	outAr = []
	curPos = 1
	maxLen = len( inAr )
	while curPos < (maxLen-binSize):
		# get subsection
		subAr = inAr[curPos:curPos+binSize+1]
		subStr = ''.join(subAr)
		count, dens = computeDensity( subStr, search )
		outAr += [ (chrmName, curPos, count, dens) ]
		curPos += binSize
	# last bin
	subAr = inAr[curPos:maxLen]
	lastC, lastD = computeDensity( subAr, search )
	outAr += [ (chrmName, curPos, lastC, lastD ) ]
	return outAr

def computeDensity( inStr, searchAr ):	
	
	lenSearch = len( searchAr[0] )
	curCount = 0
	for i in range( len(inStr) - lenSearch + 1 ):
		testStr = inStr[i:i+lenSearch]
		if testStr in searchAr:
			curCount += 1
	# end for
	testLen = float( len(inStr)-lenSearch+1 )
	# density = curCount / totalPossible
	density = (0 if testLen == 0 else float( curCount ) / testLen )
	return curCount, density

def writeOutput( outFileStr, outAr ):
	outFile = open( outFileStr, 'w' )
	outFile.write( 'chrm\tpos\tcount\tdensity\n' )
	for tup in outAr:
		outFile.write( '{:s}\t{:d}\t{:d}\t{:f}\n'.format( tup[0], tup[1], tup[2],tup[3] ) )
	outFile.close()

def parseInputs( argv ):
	binSize = BINSIZE
	outID = None
	startInd = 0
	
	for i in range(min(2,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-b=' ):
			try:
				binSize = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: bin size must be integer' )
				exit()
	# end for
	nSearch = argv[startInd]
	fastaFileStr = argv[startInd+1]
	
	if outID == None:
		outID = bth_util.fileBaseName( fastaFileStr )
	
	processInputs( fastaFileStr, nSearch, binSize, outID )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3 fasta_nuc_density [-b=bin_size] [-o=outID] <nuc_search> <fasta_file>")
	else:
		parseInputs( sys.argv[1:] )
