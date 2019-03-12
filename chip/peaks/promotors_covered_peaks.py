import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 promotors_covered_peaks.py [-gene] [-f=fasta_index] [-m=faction_covered] [-di=min_promotor_dist] [-dx=max_promotor_dist] [-o=out_id]  <gff_file> <peak_bed_file>

FRAC=2/3
MINDIST=100
MAXDIST=1000
NBINS=20

def processInputs( gffFileStr, peakBedStr, fastaIndexStr, outID, fractionCovered, includeGene, minDist, maxDist ):
	print( 'Reading GFF...' )
	gffDict, chrmDict = readGFF( gffFileStr, includeGene )
	
	if fastaIndexStr != None:
		chrmDict = readFastaIndex( fastaIndexStr )
	# prep peak dict
	emptyDict = prepDict( chrmDict )
	chrmStr = [ '{:s}({:d})'.format( chrm, chrmDict[chrm] ) for chrm in sorted(chrmDict.keys()) ]
	print( '\nGFF File: {:s}\nPeaks File: {:s}\nChromosomes (len): {:s}\n'.format( os.path.basename( gffFileStr ), os.path.basename( peakBedStr ), ', '.join( chrmStr ) ) )
	# read peak file
	print( 'Reading peak file...' )
	peakDict = readPeakFile( peakBedStr, emptyDict )
	
	print( 'Analyzing promotors and peaks...')
	covDict = handleGenePeaks( gffDict, peakDict, fractionCovered, minDist, maxDist, includeGene )
	
	outFileStr = 'promotor_peaks_'
	if outID == None:
		bName = os.path.basename( peakBedStr )
		rInd = bName.rfind( '.' )
		outFileStr += ( bName if rInd == -1 else bName[:rInd] )
	else:
		outFileStr += outID
	outFileStr += '.tsv'
	info = "#from_script:promotors_covered_peaks.py; promotor_min_dist:{:d}; promotor_max_dist:{:d}; min_coverage:{:.4f}; peak_file:{:s}; gff_file:{:s}".format( minDist, maxDist, fractionCovered, os.path.basename( peakBedStr ), os.path.basename( gffFileStr ) )
	print( 'writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, covDict, includeGene, info )
	

def readGFF( gffFileStr, includeGene ):
	gffDict = {}
	chrmDict = {}
	currChrm = None
	currMax = 0
	
	gffFile = open( gffFileStr, 'r' )
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		if lineAr[2] != 'gene':
			continue
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		strand = lineAr[6]
		gName = searchNotes( lineAr[8], 'Name' )
		gDesc = searchNotes( lineAr[8], 'Description' )
		# new chromosome
		if currChrm != None and currChrm != chrm:
			chrmDict[ chrm ] = currMax
			currChrm = chrm
			currMax = 0
		# check max
		currMax = ( end if end > currMax else currMax )
		# add gene
		if gffDict.get( chrm ) == None:
			gffDict[ chrm ] = []
		
		if includeGene:
			gffDict[ chrm ] += [ ( start, end, gName, strand, gDesc ) ]
		else:
			gffDict[ chrm ] += [ (start, end, gName, strand) ]
	# end for line
	# save last
	if currChrm != None and currMax != 0:
		chrmDict[ currChrm ] = currMax
	gffFile.close()
	
	return gffDict, chrmDict

def searchNotes( notesStr, searchStr ):
	search = searchStr + "="
	index = notesStr.find(search)
	if index == -1:
		return ""
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def readFastaIndex( fastaIndexStr ):
	chrmDict = {}
	fastaIndex = open( fastaIndexStr, 'r' )
	for line in fastaIndex:
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) length ...
		chrmDict[ lineAr[0] ] = int( lineAr[1] )
	fastaIndex.close()
	return chrmDict

def prepDict( chrmDict ):
	outDict = {}
	for chrm in chrmDict.keys():
		outDict[ chrm ] = [0] * (chrmDict[chrm]+1)
	return outDict

def readPeakFile( peakBedStr, peakDict ):
	chrms = peakDict.keys()
	peakFile = open( peakBedStr, 'r' )
	for line in peakFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		chrm = lineAr[0]
		if chrm not in chrms:
			continue
		start = int( lineAr[1] ) + 1
		end = int( lineAr[2] ) + 1
		for i in range( start, end + 1):
			peakDict[chrm][i] = 1
	# end for line
	peakFile.close()
	return peakDict

def handleGenePeaks( gffDict, peakDict, fractionCovered, minDist, maxDist, includeGene ):
	outDict = {}
	# loop through chromosomes
	for chrm in peakDict:
		outDict[chrm] = []
		peakAr = peakDict[ chrm ]
		# loop through genes
		for gene in gffDict[ chrm ]:
			start = gene[0]
			end = gene[1]
			name = gene[2]
			strand = gene[3]
			desc = gene[4]
			if strand == "+":	# positive
				pStart = start - maxDist
				pEnd = start - minDist
			else:	# negative
				pStart = end + minDist
				pEnd = end + maxDist
			if pStart < 0 or pEnd > len( peakAr ):
				continue
			# handle promotor region
			resB, resE, resC = handlePromotorPeaks( peakAr, pStart, pEnd, fractionCovered )
			# covered promotor but don't include gene
			if resB == None or resE == None:
				continue
			elif includeGene == False:
				outDict[chrm] += [ (resB, resE, resC, strand, name) ]
			# covered promotor -> check gene
			else:
				frac = computeFractionCovered( start, end, peakAr )
				isCov = (frac >= fractionCovered)
				outDict[chrm] += [ (resB, resE, resC, name, start, end, strand, isCov, desc ) ]
		# end for gene
	# end for chrm
	return outDict

def handlePromotorPeaks( peakAr, start, end, fractionCovered ):
	
	binWidth = int( math.ceil( (end - start + 1)/NBINS) )
	beginning = None
	ending = None
	covered = 0
	total = 0
	
	# loop through bins
	for bin in range(NBINS):
		# new start and end
		tStart = start + binWidth * bin
		tEnd = min( end, tStart + binWidth )
		#print( tStart, tEnd )
		cov, tot = computeCoverage( tStart, tEnd, peakAr )
		coverage = float( cov ) / tot
		if coverage >= fractionCovered:
			# add to counts
			covered += cov
			total += tot
			# first covered region
			if beginning == None:
				beginning = tStart
				ending = tEnd
			# extend
			else:
				ending = tEnd
	# end for bin
	return beginning, ending, ( 0 if total == 0 else float(covered) / total )
		
def computeFractionCovered( start, end, chrmAr ):
	total = 0
	covered = 0
	for i in range( start, end+1 ):
		total += 1
		covered += chrmAr[i]
	return float( covered ) / total

def computeCoverage( start, end, chrmAr ):
	total = 0
	covered = 0
	for i in range( start, end ):
		total += 1
		covered += chrmAr[i]
	#print (covered, total )
	return covered, total

def writeOutput( outFileStr, covDict, includeGene, info ):
	# not includeGene: [ (resB, resE, resC, strand, name) ]
	# includeGene: [ (resB, resE, resC, name, start, end, strand, isCov, desc ) ]
	if includeGene:
		headerAr = [ 'chrm', 'promotor.start', 'promotor.end', 'promotor.cov', 'gene.name', 'gene.covered','gene.start', 'gene.end', 'strand', 'gene.desc' ]
	else:
		headerAr = [ 'chrm', 'promotor.start', 'promotor.end', 'promotor.cov', 'strand', 'gene' ]
	
	header = '\t'.join( headerAr )
	outFile = open( outFileStr, 'w' )
	outFile.write( info + "\n" + header + "\n" )
	
	for chrm in sorted( covDict.keys() ):
		# loop through covered promotors
		for i in range(len(covDict[ chrm ])):
			if includeGene:
				pStart, pEnd, pCov, gName, gStart, gEnd, strand, isCov, gDesc = covDict[chrm][i]
				#print( pStart, pEnd, pCov, gName, gStart, gEnd, strand, isCov, gDesc )
				outStr = "{:s}\t{:d}\t{:d}\t{:.4f}\t{:s}\t{:s}\t{:d}\t{:d}\t{:s}\t{:s}\n".format( chrm, pStart, pEnd, pCov, gName, str(isCov), gStart, gEnd, strand, gDesc )
				outFile.write( outStr )
			else:
				pStart, pEnd, pCov, strand, gName = covDict[chrm][i]
				outFile.write( '{:s}\t{:d}\t{:d}\t{:.4f}\t{:s}\t{:s}\n'.format( chrm, pStart, pEnd, pCov, strand, gName ) )
	# end for chrm
	outFile.close()
	
def parseInputs( argv ):
	# promotors_covered_peaks.py [-gene | -g] [-f=fasta_index] [-m=faction_covered] [-di=min_promotor_dist] [-dx=max_promotor_dist] [-o=out_id]  <gff_file> <peak_bed_file>
	
	fractionCovered = float(FRAC)
	includeGene = False
	minDist = MINDIST
	maxDist = MAXDIST
	outID = None
	fastaIndexStr = None
	startInd = 0
	
	for i in range(min(6, len(argv)-2 )):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			try:
				fractionCovered = float( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: fraction of gene covered must be numeric' )
				exit()
		elif argv[i].startswith( '-f=' ):
			fastaIndexStr = argv[i][3:]
			startInd += 1
		elif argv[i] == '-gene' or argv[i] == '-g':
			includeGene = True
			startInd += 1
		elif argv[i].startswith( '-di=' ):
			try:
				minDist = int( argv[i][4:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum promotor distance must be an integer' )
				exit()
		elif argv[i].startswith( '-dx=' ):
			try:
				maxDist = int( argv[i][4:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: maximum promotor distance must be an integer' )
				exit()
		elif argv[i].startswith( '-' ):
			fInd = argv[i].find( '=' )
			if fInd != -1:
				err = argv[i][:fInd] 
			else:
				fInd = argv[i][1:].find( '-')
				if fInd == -1:
					err = argv[i]
				else:
					err = argv[i][:(fInd+1)]
			print( 'ERROR: {:s} is not a valid option'.format(err) )
			printHelp()
			exit()
	# end for
	gffFileStr = argv[startInd]
	peakBedStr = argv[startInd+1]
	processInputs( gffFileStr, peakBedStr, fastaIndexStr, outID, fractionCovered, includeGene, minDist, maxDist )
	
def printHelp():
	# Usage: python3.4 promotors_covered_peaks.py [-gene] [-f=fasta_index] [-m=fraction_covered] [-di=min_promotor_dist] [-dx=max_promotor_dist] [-o=out_id]  <gff_file> <peak_bed_file>
	
	outStr = "Usage: python3.4 promotors_covered_peaks.py [-gene] [-f=fasta_index] [-m=faction_covered] [-di=min_promotor_dist] [-dx=max_promotor_dist] [-o=out_id]  <gff_file> <peak_bed_file>\n"
	outStr += "Required:\n" + "gff_file\tGFF formatted file containing genes\n" + "peak_bed_file\tBED formatted file containing ChIP peaks\n"
	outStr += "Optional:\n"
	outStr += "-gene\tinclude info about gene coverage in output [default: off]\n"
	outStr += "-f=fasta_index\tfasta index file of chromosomes [default: use GFF coordinates]\n"
	outStr += '-m=fraction_covered\tfraction of promotor/gene to be covered by peak [default: 2/3]\n'
	outStr += '-di=min_promotor_dist\tdistance in bp from end of promotor to TSS [default: 100]\n'
	outStr += '-dx=max_promotor_dist\tdistance in bp from start of promotors to TSS [default: 1000]\n'
	outStr += '*note: uses 20 windows to determine covered promotor region, if any\n'
	outStr += '-out_id\tid to be used for output file [default: out]\n'
	print( outStr )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
