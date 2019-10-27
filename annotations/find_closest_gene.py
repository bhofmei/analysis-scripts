import sys, math, glob, multiprocessing, subprocess, os, bisect, random, gzip
# import pysam, pybedtools

MAXDIST=2000
MINDIST=500

# Usage: python find_closest_gene.py [-h] [-q] [-x=max_dist] [-n=min_dist] [-o=out_id] <gff_file> <bed_file>

def processInputs( gffFileStr, bedFileStr, maxDist, minDist, outId, isPrint ):
	if isPrint:
		print( 'GFF file:', os.path.basename(gffFileStr) )
		print( 'BED file:', os.path.basename(bedFileStr) )
		print( 'Min distance:', minDist )
		print( 'Max distance:', maxDist )
	
	if outId == None:
		outFileStr = bedFileStr.replace('.bed', '')
		outFileStr += '_close-genes.tsv'
	else:
		outFileStr = outId + '_close-genes.tsv'
	
	if isPrint:
		print( 'Reading GFF' )
	
	gffDict = readGff( gffFileStr )
	
	if isPrint( 'Reading BED' )
	
	bedDict = readBed( bedFileStr )
	
	

def readGff( gffFileStr ):
	gffFile = open( gffFileStr, 'r' )
	
	outDict = []
	
	for line in gffFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or len(lineAr) < 8:
			continue
		chrm = lineAr[0]
		start = int(lineAr[3]) - 1
		end = int(lineAr[4])
		strand = lineAr[6]
		fType = lineAr[2]
		if fType != 'gene':
			continue
		gName = searchId (lineAr[8] )
		if outDict.get(chrm) == None:
			outDict[chrm] = []
		# add gene start
		bisect.insort( outDict[chrm], (start, 'start', gName, strand) )
		bisect.insort( outDict[chrm], (end, 'end', gName, strand) )
	# end for line
	
	gffFile.close()
	return outDict
		
def searchId( notesStr ):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex != -1:
		return notesStr[adIndex:endIndex+adIndex]
	else:
		return notesStr[adIndex:]

def readBed( bedFileStr ):
	
	bedFile = open( bedFileStr, 'r' )
	outDict = {}
	rCount = 1
	
	for line in bedFile:
		lineAr = line.rstrip().split('\t')
		if line.startswith('#') or len(lineAr) < 3:
			continue
		chrm = lineAr[0]
		start = int(lineAr[1])
		end = int( lineAr[2] )
		rName = 'region-'+str(rCount) if len(lineAr) == 3 else lineAr[3]
		rCount += 1
		if outDict.get(chrm) == None:
			outDict[chrm] = []
		outDict[chrm] += [(start, end, rName)]
	# end for line
	bedFile.close()
	return outDict
	
def parseInputs( argv ):
	maxDist = MAXDIST
	minDist = MINDIST
	outId = None
	isPrint = True
	startInd += 1
	
	for i in range(len(argv)):
		if argv[i] in [ '-h', '-help', '--help' ]:
			printHelp()
			exit()
		elif argv[i] == '-q':
			isPrint = True
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outId = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-x=' ):
			try:
				maxStr = (argv[i][3:]).lower()
				if maxStr.endswith( 'k' ) or maxStr.endswith( 'kb' ):
					maxDist = int( maxStr[:-1] ) * 1000
				elif maxStr.endswith( 'm' ) or maxStr.endswith( 'mb' ):
					maxDist = int( maxStr[:-1] ) * 1000000
				else:
					maxDist = int( maxStr )
				
			except ValueError:
				print( 'WARNING: max distance must be an integer...using default', MAXDIST )
				maxDist = MAXDIST
			startInd += 1
		elif argv[i].startswith( '-n=' ):
			try:
				minStr = (argv[i][3:]).lower()
				if minStr.endswith( 'k' ):
					minDist = int( minStr[:-1] ) * 1000
				elif minStr.endswith( 'm' ):
					minDist = int( minStr[:-1] ) * 1000000
				else:
					minDist = int( minStr )
				
			except ValueError:
				print( 'WARNING: min distance must be an integer...using default', MINDIST )
				minDist = MINDIST
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for i
		
	gffFileStr = argv[startInd]
	bedFileStr = argv[startInd + 1]
	
	processInputs( gffFileStr, bedFileStr, maxDist, minDist, outId, isPrint )

def printHelp():
	print('Usage: python find_closest_gene.py [-h] [-q] [-x=max_dist] [-n=min_dist] <gff_file> <bed_file>')
	
if __name__ == "__main__":
	if len(sys.argv) < 3:
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
