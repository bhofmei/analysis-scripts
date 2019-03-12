import sys, math, glob, multiprocessing, subprocess, os, bisect

# Usage: analyze_smp_location.py [-o=out_prefix] <gff_file> <smp_file> <gBM_file>
# Using the GFF file for annotation information and gBM file to separate gBM and non-gBM genes, outputs number of smps for each category: gBM genes CDS, gBM genes introns, gBM genes UTRs, non-GBM genes CDS, non-GBM genes introns, non-GBM genes UTRs, TEs, other
# note: only using first transcript of each gene

def processInputs( gffFileStr, smpFileStr, outFileStr, outCountStr ):
	print( 'Reading SMP file...' )
	smpDict, smpCount, header = readSMP( smpFileStr )
	print( 'Reading GFF file...' )
	featureAr = readGFF( gffFileStr )
	print( 'Sorting and counting SMPs...' )
	countAr = countSMPs( smpDict, featureAr )
	print( 'Writing output to {:s}...'.format( outCountStr ) )
	writeOutput( outCountStr, countAr, smpCount )

def readSMP( smpFileStr ):
	
	smpFile = open( smpFileStr, 'r' )
	smpDict = {}
	# key is chromosome and value is bisect array of positions
	smpCount = 0
	header = ''
	
	for line in smpFile:
		if line.startswith( '#' ):
			header += line
			continue
		elif line == '\n':
			continue
		lineAr = line.rstrip().split( '\t' )
		
		chrm = lineAr[0]
		pos = int( lineAr[1] )
		# chrm not in dict
		if smpDict.get( chrm ) == None:
			smpDict[ chrm ] = []
		bisect.insort( smpDict[ chrm ], pos )
		smpCount += 1
	smpFile.close()
	return smpDict, smpCount, header

def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

def readGFF( gffFileStr, gbmAr ):
	
	gffFile = open( gffFileStr, 'r' )
	featureAr = []
	prevPos = None
	
	for line in gffFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = min(int(lineAr[3]),int(lineAr[4]))
		end = max(int(lineAr[3]),int(lineAr[4]))
		chrm = lineAr[0]
		strand = lineAr[6]
		
		# handle mRNA for intron tracking
		if lineAr[2] == 'mRNA':
			gName, isPrimary = getGene( lineAr[8], 'Name' )
			if isPrimary == False: # don't want secondary transcripts
				prevPos = None
		# reset intron tracking at gene
		elif lineAr[2] == 'gene':
			prevPos = None
		# handle TEs
		elif lineAr[2] == 'transposable_element':
			featureAr += [ (chrm, start, end, 'TE') ]
		# handle CDS and UTRs ( with introns )
		elif lineAr[2] in [ 'CDS','five_prime_UTR','three_prime_UTR'] :
			gName, isPrimary = getGene( lineAr[8], 'Parent' )
			fType = ('UTR' if 'UTR' in lineAr[2] else 'CDS' )
			if isPrimary:
				# positive strand introns
				if strand == '+':
					# no prevPos
					if prevPos == None or prevPos == start:
						prevPos = end + 1
					else:
						fTypeI = 'intron'
						featureAr += [ (chrm, prevPos, start-1, fTypeI) ]
						prevPos = end + 1
				elif strand == '-':
					if prevPos == None or prevPos == end:
						prevPos = start - 1
					else:
						fTypeI = 'intron'
						featureAr += [ (chrm, end+1, prevPos, fTypeI) ]
						prevPos = start - 1
				featureAr += [ (chrm, start, end, fType) ]
	# end for
	gffFile.close()
	return featureAr

def getGene( notesStr, searchStr ):
	search = searchStr + '='
	index = notesStr.find(search)
	adIndex = index + len(search)
	#try comma
	cIndex = notesStr[adIndex:].find(',')
	sIndex = notesStr[adIndex:].find(';')
	if cIndex == -1 and sIndex == -1:
		name = notesStr[ adIndex: ]
	elif cIndex == -1:
		name = notesStr[adIndex:adIndex+sIndex]
	else:
		name = notesStr[adIndex:adIndex+cIndex]
	isPrimary = name.endswith( '.1' )
	pInd = name.rfind('.')
	return name[:pInd], isPrimary

def countSMPs( smpDict, featureAr ):
	
	fTypeAr = [ 'CDS','UTR','intron','TE']
	countAr = [0] * (len(fTypeAr))
	# (0) CDS (1) UTR (2) intron (3) TE (4) other
	
	# loop through genes
	for feat in featureAr:
		chrm = feat[0]
		start = feat[1]
		end = feat[2]
		fType = feat[3]
		try:
			smpAr = smpDict[ chrm ]
			fInd = fTypeAr.index( fType )
		except KeyError:
			continue
		except ValueError:
			fInd = len( countAr )
		for i in range( start, end+1 ):
			isIn = bisectIndex( smpAr, i )
			if isIn != None:
				countAr[fInd] += 1
				
	return countAr

def writeOutput( outFileStr, countAr, smpCount ):
	
	outFile = open( outFileStr, 'w' )
	fTypeAr = ['CDS', 'UTR', 'intron', 'TE' ,'other' ]
	updateCountAr = countAr + [ smpCount - sum(countAr) ]
	for i in range(len(updateCountAr)):
		outFile.write( '{:s}\t{:d}\n'.format( fTypeAr[i], updateCountAr[i] ) )
	outFile.write( 'total\t{:d}\n'.format( smpCount ) )
	outFile.close()

def parseInputs( argv ):
	
	startInd = 0
	outPre = None
	
	if argv[0].startswith( 'o=' ):
		outPre = argv[0][3:]
		startInd += 1
	
	gffFileStr = argv[startInd]
	smpFileStr = argv[startInd+1]
	
	if outPre == None:
		baseName = os.path.basename(smpFileStr)
		outPre = baseName[:(baseName.rfind('.'))]
	outCountStr = outPre + '_loc_count.tsv'
	outFileStr = outPre + '_loc.tsv'
	processInputs( gffFileStr, smpFileStr, outFileStr, outCountStr )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: analyze_smp_location.py [-o=out_prefix] <gff_file> <smp_file>")
	else:
		parseInputs( sys.argv[1:] )
