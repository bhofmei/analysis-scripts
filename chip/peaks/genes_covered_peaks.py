import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 genes_covered_peaks.py [-desc] [-f=fasta_index] [-m=faction_gene_covered] [-o=out_id]  <gff_file> <peak_bed_file>

FRAC=2/3

def processInputs( gffFileStr, peakBedStr, fastaIndexStr, outID, fractionCovered, includeDesc ):
	info = '#from_script:genes_covered_peaks.py; gff_file:{:s}; peak_bed_file:{:s}'.format( os.path.basename( gffFileStr ), os.path.basename( peakBedStr ) )
	print( 'Reading GFF...' )
	gffDict, chrmDict = readGFF( gffFileStr, includeDesc )
	
	if fastaIndexStr != None:
		chrmDict = readFastaIndex( fastaIndexStr )
		info += ';fasta_index_file:{:s}'.format( os.path.basename( fastaIndexStr ) )
	# prep peak dict
	emptyDict = prepDict( chrmDict )
	chrmStr = [ '{:s}({:d})'.format( chrm, chrmDict[chrm] ) for chrm in sorted(chrmDict.keys()) ]
	print( '\nGFF File: {:s}\nPeaks File: {:s}\nChromosomes (len): {:s}\n'.format( os.path.basename( gffFileStr ), os.path.basename( peakBedStr ), ', '.join( chrmStr ) ) )
	# read peak file
	print( 'Reading peak file...' )
	peakDict = readPeakFile( peakBedStr, emptyDict )
	info += ';fraction_covered_cutoff:{:.4f}'.format( fractionCovered )
	print( 'Analyzing genes and peaks...')
	covDict = handleGenesPeaks( gffDict, peakDict, fractionCovered )
	
	outFileStr = 'peaks_'
	if outID == None:
		bName = os.path.basename( peakBedStr )
		rInd = bName.rfind( '.' )
		outFileStr += ( bName if rInd == -1 else bName[:rInd] )
	else:
		outFileStr += outID
	outFileStr += '.tsv'
	
	print( 'writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, covDict, includeDesc, info )
	

def readGFF( gffFileStr, includeDesc ):
	gffDict = {}
	chrmDict = {}
	currChrm = None
	currMax = 0
	
	gffFile = open( gffFileStr, 'r' )
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) < 3:
			print(lineAr)
			continue
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		if lineAr[2] != 'gene':
			continue
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
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
		if includeDesc:
			gffDict[ chrm ] += [ ( start, end, gName, gDesc ) ]
		else:
			gffDict[ chrm ] += [ (start, end, gName) ]
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

def handleGenesPeaks( gffDict, peakDict, fractionCovered ):
	outDict = {}
	
	# loop through chrms
	for chrm in peakDict:
		outDict[chrm] = []
		peakAr = peakDict[ chrm ]
		# loop through genes
		for gene in gffDict[chrm]:
			start = gene[0]
			end = gene[1]
			name = gene[2]
			frac = computeFractionCovered( start, end, peakAr )
			if frac >= fractionCovered:
				if len(gene) > 3:
					outDict[ chrm ] += [ (start, end, name, gene[3], frac ) ]
				else:
					outDict[ chrm ] += [ (start, end, name, frac ) ]
		# end for gene
	# end for chrm
	return outDict

def computeFractionCovered( start, end, chrmAr ):
	total = 0
	covered = 0
	for i in range( start, end+1 ):
		total += 1
		covered += chrmAr[i]
	return float( covered ) / total

def writeOutput( outFileStr, covDict, includeDesc, info ):
	if includeDesc:
		header = 'chrm\tgene.start\tgene.end\tgene.name\tgene.description\tfrac.cov\n'
	else:
		header = 'chrm\tgene.start\tgene.end\tgene.name\tfrac.cov\n'
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' + header )
	
	for chrm in sorted( covDict.keys() ):
		# loop through covered genes
		for i in range(len(covDict[ chrm ])):
			if includeDesc:
				start, end, name, desc, frac = covDict[chrm][i]
				outFile.write( '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:.4f}\n'.format( chrm, start, end, name, desc, frac ) )
			else:
				start, end, name, frac = covDict[chrm][i]
				outFile.write( '{:s}\t{:d}\t{:d}\t{:s}\t{:.4f}\n'.format( chrm, start, end, name, frac ) )
	# end for chrm
	outFile.close()
	
def parseInputs( argv ):
	fractionCovered = float(FRAC)
	outID = None
	fastaIndexStr = None
	startInd = 0
	includeDesc = False
	
	for i in range(min(4, len(argv)-2 )):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			try:
				fractionCovered = float( argv[i][3:] )
				if fractionCovered > 1:
					fractionCovered = fractionCovered / 100.0
				startInd += 1
			except ValueError:
				print( 'ERROR: fraction of gene covered must be numeric' )
				exit()
		elif argv[i].startswith( '-f=' ):
			fastaIndexStr = argv[i][3:]
			startInd += 1
		elif argv[i] == '-d':
			includeDesc = True
			startInd += 1
	# end for
	gffFileStr = argv[startInd]
	peakBedStr = argv[startInd+1]
	processInputs( gffFileStr, peakBedStr, fastaIndexStr, outID, fractionCovered, includeDesc )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 genes_covered_peaks.py [-m=faction_gene_covered] [-o=out_id] [-f=fasta_index] <gff_file> <peak_bed_file>")
	else:
		parseInputs( sys.argv[1:] )
