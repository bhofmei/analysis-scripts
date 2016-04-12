import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 genes_in_range_bed.py <bed_file> <gff_file> <distance> [out_prefix]

def readGFF(gffFileStr):
	
	gffFile = open (gffFileStr, 'r' )
	geneDict = {}
	count = 0
	for line in gffFile:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		lineAr = line.split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		scaffold = lineAr[0]
		strand = lineAr[6]
		#print ( "GFF:",lineAr)
		
		# only apply to type = gene
		if lineAr[2] == "gene":
			count += 1
			name = getGeneName( lineAr[8] )
			# put into geneDict
			if geneDict.get(scaffold) == None:
				geneDict[scaffold] = []
			geneDict[scaffold] += [(start, end, name)]
	
	gffFile.close()
	#print( count )
	return geneDict, count

def getGeneName (notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def readBed( bedFileStr, distance ):
	# note these aren't actually full bed files
	bedFile = open( bedFileStr, 'r' )
	bedDict = {}
	
	for line in bedFile:
		lineAr = line.rstrip().split()
		# (0) chrm (1) start (2) end
		chrm = lineAr[0]
		start = int( lineAr[1] )
		end = int( lineAr[2] )
		if bedDict.get( chrm ) == None:
			bedDict[chrm] = []
		bedDict[chrm] += [ (start-distance,end+distance) ]
	bedFile.close()
	return bedDict

def determineOverlap( bedDict, gffDict ):
	# bed dict: {chrm:[(start,end),(start,end),...],chrm2:...}
	# gene dict: {chrm:[(start,end,name),(start,end,name)],chrm2:...}

	overlapAr = []
	for chrm in gffDict.keys():
		if bedDict.get( chrm ) == None:
			continue
		regionAr = sorted( bedDict.get( chrm ) )
		geneAr = sorted( gffDict.get( chrm ) )
		currPos = 0
		currRegion, rStart, rEnd = getRegion( regionAr, currPos )
		# loop through genes
		for i in range( len( geneAr ) ):
			gStart = geneAr[i][0]
			gEnd = geneAr[i][1]
			#print( currRegion, geneAr[i] )
			while rEnd < gStart and currPos < len( regionAr )-1:
				currPos += 1
				currRegion, rStart, rEnd = getRegion( regionAr, currPos )
			# gene starts before region
			if gStart < rStart and rStart < gEnd:
				overlapAr += [ geneAr[i][2] ]
				#print( 'True' )
			# gene entirely within region
			elif rStart <= gStart and gEnd <= rEnd:
				overlapAr += [ geneAr[i][2] ]
				#print( 'True' )
			# gene ends after region
			elif gStart < rEnd and rEnd < gEnd:
				overlapAr += [ geneAr[i][2] ]
				#print( 'True' )
		# finish looping through genes on chrm
	# finished will all chromosomes
	return overlapAr

def getRegion( regionAr, pos ):
	region = regionAr[pos]
	return region, region[0], region[1]

def writeOutput( outFileStr, overlapAr, headerStr ):
	outFile = open( outFileStr, 'w' )
	outFile.write( headerStr + '\n' )
	outFile.write( '\n'.join( overlapAr ) + '\n' )
	outFile.close()

def processInputs( bedFileStr, gffFileStr, distance, outPre ):
	outFileStr = outPre + '_' + str( distance ) + '.txt'
	print( 'Reading', bedFileStr, '...' )
	bedDict = readBed( bedFileStr, distance*1000 )
	#print( bedDict )
	print( 'Reading', gffFileStr, '...' )
	gffDict, geneCount = readGFF( gffFileStr )
	#print( gffDict )
	print( 'Determining overlap...' )
	overlapAr = determineOverlap( bedDict, gffDict )
	headerStr = '# genes within {:d} kb of region in {:s} - {:d} of {:d} genes'.format( distance, bedFileStr, len( overlapAr ), geneCount )
	print( 'Writing output to', outFileStr, '...' )
	writeOutput( outFileStr, overlapAr, headerStr )
	print( 'Done.' )

def parseInputs( argv ):
	bedFileStr = argv[0]
	gffFileStr = argv[1]
	try:
		distance = int( argv[2] )
	except ValueError:
		print( 'ERROR: distance must be integer.' )
		exit()
	if len( argv ) == 4:
		outPre = argv[3]
	else:
		rind = bedFileStr.rindex( '.' )
		outPre = bedFileStr[:rind] + '_genes'
	processInputs( bedFileStr, gffFileStr, distance, outPre ) 

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		print ("Usage: python3.4 genes_in_range_bed.py <bed_file> <gff_file> <distance> [out_prefix]")
	else:
		parseInputs( sys.argv[1:] )
