import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: python3 portion_from_bed_pe.py [-genes|-cds] [-p=num_proc] [-ex=exclude_genes|-in=include_genes] [-f=fasta_index] <gff_file> <bed_file> [bed_file]
# computes portion of reads in BED file that overlap with a particular feature type

FEATURES='gene'
NUMPROC=1

def processInputs( bedFileAr, gffFileStr, featureType, numProc, excludeList, includeList, fastaIndex ):
	
	info = '#from_script:portion_from_bed_pe.py; gff_file: {:s}; features: {:s}; '.format( os.path.basename( gffFileStr ), featureType )
	if includeList != None:
		info += 'include_file: {:s}; '.format( includeFile )
	elif excludeList != None:
		info += 'exclude_file: {:s}; '.format( excludeFile )
	smFileAr = [os.path.basename(x) for x in bedFileAr ]
	info += 'bed_files: {:s}'.format( ','.join( smFileAr ) )
	
	print( 'GFF File: {:s}\nBED Files: {:s}\nInclude File: {:s}\nExclude File: {:s}\nFeatures: {:s}'.format( os.path.basename( gffFileStr ), ','.join( smFileAr ), str( includeList ), str( excludeList ), featureType ) )

	# read include/exclude files
	if includeList != None:
		includeFile = includeList
		includeList = readBasicFile( includeFile )
	elif excludeList != None:
		excludeFile = excludeList
		excludeList = readBasicFile( excludeFile )
	
	# read fasta index if it exists
	if fastaIndex == None:
		chrmDict = None
	else:
		print( 'Reading FASTA index' )
		fiFile = FileFASTAIndex( fastaIndex )
		chrmDict = fiFile.getLengths()
	
	# read gff with include/exclude information and feature type
	print( 'Reading GFF' )
	gffDict, nfeat = readGFF( gffFileStr, chrmDict, featureType, excludeList, includeList )
	print( 'Number of features in GFF: {:d}'.format( nfeat ) )

	# read BED files and process
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( readBED, args=(f, gffDict ) ) for f in bedFileAr ]
	resultsAr = [ p.get() for p in results ]
	
	# write output
	outFileStr = 'out_portion_{:s}.tsv'.format( featureType )
	print( 'Writing output to {:s}'.format( outFileStr ) )
	writeOutput( resultsAr, bedFileAr, outFileStr, info )
	print( 'Done' )
	
def readBasicFile( fileStr ):
	# this file must use # to indicate comments
	# for this function, there is only one column of information
	print( 'Reading file {:s}'.format( os.path.basename( fileStr ) ) )
	inFile = open( fileStr, 'r' )
	outAr = []
	for line in inFile:
		l = line.rstrip().lstrip()
		if l.startswith( '#' ) == False:
			outAr += [ l ]
	inFile.close()
	return outAr

def readGFF( gffFileStr, chrmDict, featureType, excludeList, includeList ):
	if chrmDict == None:
		outDict, n = readGFFNoLenths( gffFileStr, featureType, excludeList, includeList )
	else:
		outDict, n = readGFFLengths( gffFileStr, chrmDict, featureType, excludeList, includeList )
	
	return outDict, n

def readGFFNoLenths( gffFileStr, featureType, excludeList, includeList ):
	print( 'ERROR: coming soon...for now specify the fasta index with the -f option' )
	exit()

def readGFFLengths( gffFileStr, chrmDict, featureType, excludeList, includeList ):
	outDict = {}
	count = 0
	# prepare dictionary
	for chrm in chrmDict.keys():
		cLen = chrmDict[chrm]
		outDict[chrm] = [False] * cLen
	
	# read GFF file
	gffFile = open( gffFileStr, 'r' )
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		chrm = lineAr[0]
		# chrm doesn't exist or not correct feature
		if chrmDict.get( chrm ) == None or lineAr[2] != featureType:
			continue
			
		start = int( lineAr[3] ) - 1
		end = int( lineAr[4] )
		gName = getGeneName( lineAr[8] )
		
		# not part of include list
		if includeList != None and gName not in includeList:
			continue
		# is part of exclude list
		elif excludeList != None and gName in excludeList:
			continue
		# otherwise keep it	
		else:
			count += 1
			try:
				for i in range(start, end ):
					outDict[chrm][i] = True 
			except IndexError:
				print( 'WARNING: feature {:s} goes farther than chromosome is long ({:d} vs {:d})'.format( gName, end, len(outDict[chrm] ) ) )	
	# end for line
	return outDict, count

def getGeneName( notesStr ):
	# if it has Parent= use that (CDS) as gene name otherwise use Name= (gene)
	search1 = "Parent="
	search2 = "Name="
	index = notesStr.find(search1)
	adIndex = index + len(search1)
	if index == -1:
		index = notesStr.find( search2 )
		adIndex = index + len( search2 )
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return notesStr[adIndex:endIndex+adIndex]

def readBED( bedFileStr, gffDict ):
	print( 'Reading {:s}'.format( os.path.basename( bedFileStr ) ) )
	resultsDict = {}
	bedFile = open( bedFileStr, 'r' )
	tmpTrue = 0
	tmpTotal = 0
	curChrm = None
	curChrmAr = []
	totalTrue = 0
	totalTotal = 0
	print( list(gffDict.keys() ))
	
	for line in bedFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		chrm = lineAr[0]
		start = int( lineAr[1] )-1
		end = int( lineAr[2] )
		if chrm not in gffDict.keys():
			continue
		# reset tmp chrm
		if chrm != curChrm:
			if curChrm != None:
				resultsDict[ curChrm ] = ( tmpTrue, tmpTotal )
				totalTrue += tmpTrue
				totalTotal += tmpTotal
			curChrm, tmpTrue, tmpTotal = chrm, 0, 0
			print( 'change chrm', curChrm )
			curChrmAr = gffDict[ chrm ]
		# process line
		# if it has any overlap with feature we will include it
		try:
			tmpTotal += 1
			for i in range(start, end ):
				if curChrmAr[i]:
					tmpTrue += 1
					break
		except IndexError:
			print( 'WARNING: read out of index' )
	# end for line
	bedFile.close()
	
	# clean up results
	resultsDict[ curChrm ] = ( tmpTrue, tmpTotal )
	totalTrue += tmpTrue
	totalTotal += tmpTotal
	resultsDict[ 'total' ] = ( totalTrue, totalTotal )
	return resultsDict

def writeOutput( resultsAr, bedFileAr, outFileStr, info ):
	outFile = open( outFileStr, 'w' )
	headerAr = ['sample','chrm','reads_overlap','reads_total','percentage']
	
	outFile.write(info + '\n' + '\t'.join(headerAr) + '\n')
	# loop through files
	for i in range(len(bedFileAr)):
		smName = os.path.basename(bedFileAr[i]).replace('.bed','')
		resultsDict = resultsAr[i]
		# loop through chrms
		for chrm in sorted( resultsDict.keys() ):
			tmpTrue = resultsDict[chrm][0]
			tmpTotal = resultsDict[chrm][1]
			tmpP = ('NA' if tmpTotal == 0 else '{:.2f}'.format(float(tmpTrue) / float(tmpTotal)*100))
			outStr = '{:s}\t{:s}\t{:d}\t{:d}\t{:s}\n'.format( smName, chrm, tmpTrue, tmpTotal, tmpP )
			outFile.write( outStr )
		# end for chrm
	# end for i
	outFile.close()
		
def parseInputs( argv ):
	featureType = None
	numProc = NUMPROC
	excludeList = None
	includeList = None
	fastaIndex = None
	startInd = 0
	
	for i in range(min(6,len(argv))):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif (argv[i] == '-genes' or argv[i] == '-gene'):
			if featureType != None:
				print( 'ERROR: do not specify more than one feature type' )
				exit()
			else:
				featureType = 'gene'
				startInd += 1
		elif (argv[i] == '-cds' or argv[i] == '-cd'):
			if featureType != None:
				print( 'ERROR: do not specify more than one feature type' )
				exit()
			else:
				featureType = 'CDS'
				startInd += 1
		elif argv[i].startswith( '-ex=' ):
			if includeList != None:
				print( 'ERROR: do not specify include and exclude list' )
				exit()
			else:
				excludeList = argv[i][4:]
				startInd += 1
		elif argv[i].startswith( '-in' ):
			if excludeList != None:
				print( 'ERROR: do not specify include and exclude list' )
				exit()
			else:
				includeList = argv[i][4:]
				startInd += 1
		elif argv[i].startswith( '-f=' ):
			fastaIndex = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-' ):
			print( 'WARNING: ignoring parameter {:s}'.format( argv[i] ) )
	# end for argv
	if featureType == None:
		featureType = FEATURES
	
	gffFileStr = argv[startInd]
	bedFileAr = []
	for j in range(startInd+1, len(argv)):
		bedFileAr += [ argv[j] ]
	
	processInputs( bedFileAr, gffFileStr, featureType, numProc, excludeList, includeList, fastaIndex )

def printHelp():
	print( 'Usage: python3 portion_from_bed _pe.py [-genes|-cds] [-p=num_pric] [-ex=exclude_genes|-in=include_genes] [-f=fasta_index] <gff_file> <bed_file> [bed_file]' )
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
