import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *
import bth_util

# Usage: python rpm_by_region.py [-strand | -cssr | -cssr=distance] [-o=out_id] [-p=num_proc] <region_file> <bed_file> [bed_file]*

NUMPROC=1
CSSRDIST=10000

def processInputs( regionFileStr, bedFileAr, numProc, isStrand, isCSSR, cssrDist, outID ):
	sampleNamesAr = getSampleNames( bedFileAr )
	# read region file
	print( 'Reading region file {:s}'.format( os.path.basename( regionFileStr ) ) )	
	if isCSSR:
		regionAr = readCSSRFile( regionFileStr, cssrDist )
	else:
		regionAr = readRegionFile( regionFileStr, isStrand )
	
	# process BED files
	useStrand = isStrand or isCSSR
	if outID == None:
		outID = bth_util.fileBaseName( regionFileStr )
	
	print( 'Begin processing with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processBedFile, args=(f, regionAr, useStrand) ) for f in bedFileAr ]
	outDictMat = [ p.get() for p in results ]
	
	if isCSSR:
		outFileStr = '{:s}_rpm_cssr_{:s}.tsv'.format( outID, bth_util.binSizeToStr( cssrDist ) )
	elif isStrand:
		outFileStr = '{:s}_rpm_stranded.tsv'.format( outID )
	else:
		outFileStr = '{:s}_rpm.tsv'.format( outID )
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, regionAr, outDictMat, sampleNamesAr, useStrand )
	print( 'Done' )

def getSampleNames( fileStrAr ):
	sampleNamesAr = []
	for fileStr in fileStrAr:
		leftIndex = fileStr.rfind('/')
		rightIndex = fileStr.rfind('.')
		sampleName = fileStr[leftIndex+1:rightIndex]
		sampleNamesAr += [ sampleName ]
	return sampleNamesAr
	
def readCSSRFile( inFileStr, cssrDist ):
	outAr = []
	inFile = open( inFileStr, 'rt' )
	for line in inFile:
		if line.startswith( '#' ):
			continue # comments
		# (0) cSSR name (1) chrm (2) start (3) end
		# need to make 4 separate regions actually
		# top sense, top antisense, bottom sense, bottom antisense
		lineAr = line.rstrip().split( '\t' )
		cssrName = lineAr[0]
		chrm = lineAr[1]
		start = int( lineAr[2] )
		end = int( lineAr[3] )
		sStart = max(1, start - cssrDist)
		sEnd = end + cssrDist
		regionNames = [ cssrName+x for x in ['_top_sense', '_bottom_antisense', '_top_antisense', '_bottom_sense' ] ]
		regionStarts = [sStart, sStart, end, end ]
		regionEnds = [ start, start, sEnd, sEnd ]
		regionStrand = ['+', '-', '+', '-']
		for i in range(len(regionNames)):
			outAr += [ (regionNames[i], chrm, regionStarts[i], regionEnds[i], regionStrand[i] ) ]
	# end for line
	return outAr

def readRegionFile( inFileStr, isStrand ):
	outAr = []
	inFile = open( inFileStr, 'rt' )
	for line in inFile:
		if line.startswith( '#' ):
			continue	# comment line
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) start (2) end (3?) strand
		if len(lineAr) < 4 and isStrand:
			print( 'ERROR: This line does not specify a strand. Correct the file or run again without the strand flag.\n{:s}'.format( line.rstrip() ) )
		elif len(lineAr ) < 4:	# unstranded
			chrm = lineAr[0]
			start = int( lineAr[1] )
			end = int( lineAr[2] )
			strand = ''
			name = '{:s}:{:d}-{:d}'.format( chrm, start, end )
			outAr += [ (name, chrm, start, end, strand) ]
		else:	# stranded
			chrm = lineAr[0]
			start = int( lineAr[1] )
			end = int( lineAr[2] )
			strand = lineAr[3]
			name = '{:s}:{:d}-{:d}({:s})'.format( chrm, start, end, strand )
			outAr += [ (name, chrm, start, end, strand) ]
	# end for line
	inFile.close()
	return outAr

def processBedFile( bedFileStr, regionAr, useStrand ):
	bedFile = FileBED( bedFileStr )
	bedDict, readCounts = bedFile.getBedDict( middle = True, stranded = useStrand )
	milCounts = readCounts / 1000000
	
	outDict = {}
	# loop through regions
	for region in regionAr:
		name, chrm, start, end, strand = region
		sChrm = chrm + strand
		chrmDict = bedDict.get( sChrm )
		if chrmDict != None:
			counts = countRegion( chrmDict, start, end )
			rpm = float(counts) / (milCounts)
			outDict[name] = rpm
	# end for
	return outDict		
		
def countRegion( bedDict, start, end ):
	count = 0
	for pos in range(start, end+1):
		try:
			dictEntry = bedDict.get(pos)
			if dictEntry != None:
				count += dictEntry
		except IndexError:
			continue
	# end for pos
	return count

def writeOutput( outFileStr, regionAr, outDictMat, sampleNamesAr, useStrand):
	outFile = open( outFileStr, 'w' )
	if useStrand:
		headerAr = ['#sample', 'region_name', 'chrm', 'strand', 'start', 'end', 'rpm' ]
	else:
		headerAr = ['#sample', 'region_name', 'chrm', 'start', 'end', 'rpm' ]
	outFile.write( '\t'.join(headerAr) +'\n')
	# loop through samples
	for i in range(len(sampleNamesAr)):
		sampDict = outDictMat[i]
		# loop through regions
		for region in regionAr:
			name, chrm, start, end, strand = region
			rpm = sampDict.get(name)
			if rpm == None:
				rpmS = 'NA'
			else:
				rpmS = str(rpm)
			outStr = '{:s}\t{:s}\t{:s}\t'.format( sampleNamesAr[i], name, chrm )
			if useStrand:
				outStr += '{:s}\t'.format( strand )
			outStr += '{:d}\t{:d}\t{:s}\n'.format( start, end, rpmS )
			outFile.write( outStr )
		# end for region
	# end for sample
	outFile.close()
			
		
def parseInputs( argv ):
	numProc = NUMPROC
	outID = None
	isStrand = False
	isCSSR = False
	cssrDist = CSSRDIST
	startInd = 0
	
	for i in range(min(4,len(argv))):
		if argv[i] == '-strand':
			isStrand = True
			startInd += 1
		elif argv[i].startswith('-cssr='):
			isCSSR = True
			cssrDist = bth_util.strToDistance( argv[i][6:] )
			if cssrDist == False:
				print( 'ERROR: cSSR distance specification is incorrect' )
				exit()
			startInd += 1
		elif argv[i] == '-cssr':
			isCSSR = True
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	regionFileStr = argv[startInd]
	bedFileAr = []
	for j in range(startInd+1, len(argv)):
		bedFileAr += [ argv[j] ]
	processInputs( regionFileStr, bedFileAr, numProc, isStrand, isCSSR, cssrDist, outID )

def printHelp():
	print( "Usage: python rpm_by_region.py [-strand | -cssr | -cssr=distance] [-o=out_id] [-p=num_proc] <region_file> <bed_file> [bed_file]*" )
	print( 'Required:' )
	print( 'region_file\tchrm, start, end, optional strand; tab-delim' )
	print( 'bed_file\tBED-formatted file' )
	print( 'Optional:' )
	print( '-strand\t\tstrand-specific rpm [default off]' )
	print( '-cssr\t\tconvergent strand [default distance 10kbp]' )
	print( '-cssr=N\t\tcssr with distance N' )
	print( '-o=out_id\toutput identifier' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )
	
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
