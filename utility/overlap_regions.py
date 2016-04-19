import sys

# Usage: python3.4 overlap_regions.py <regions_file_1> <regions_file_2> <out_prefix>

def readFile1( inFileStr ):
	outDict = {}
	
	inFile = open( inFileStr, 'r' )
	
	for line in inFile:
		lineAr = line.rstrip().split()
		# (0) chrm (1) start (2) end
		chrm = lineAr[0]
		start = int( lineAr[1] )
		end = int (lineAr[2] )
		
		# chrm not in dict already
		if outDict.get( chrm ) == None:
			outDict[chrm] = [ (start, end) ]
		else:
			outDict[chrm] += [ (start, end) ]
			
	inFile.close()
	return outDict

def readFile2( inFileStr, compDict ):
	
	outDict = {}
	inFile = open( inFileStr, 'r' )
	
	for line in inFile:
		lineAr = line.rstrip().split()
		chrm = lineAr[0]
		lStart = int( lineAr[1] )
		lEnd = int( lineAr[2] )
		
		chrmAr = compDict.get( chrm )
		# chrm isn't in compDict
		if chrmAr == None:
			continue
		# loop through possible regions
		for region in chrmAr:
			rStart = region[0]
			rEnd = region[1]
			
			# too far - don't waste time
			if rStart > lEnd and rEnd > lEnd:
				continue
			# in region
			elif ( rStart <= lStart and rEnd > lStart ) or ( lStart <= rStart and lEnd > rStart ):
				# chrm not in outDict
				if outDict.get( chrm ) == None:
					outDict[chrm] = set()
				outDict[chrm].add( (rStart, rEnd ) ) # coordinates of first
	
	inFile.close()
	return outDict

def writeOutput( outFileStr, outDict ):
	
	outFile = open( outFileStr, 'w' )
	
	for chrm in sorted(outDict.keys()):
		chrmSet = outDict[chrm]
		for region in chrmSet:
			outFile.write( "{:s}\t{:d}\t{:d}\n".format( chrm, region[0], region[1] ) )
			
	outFile.close()

def mainFunction( fileStr1, fileStr2, outPre):
	
	compDict = readFile1( fileStr1 )
	print( "read", fileStr1 )
	outDict = readFile2( fileStr2, compDict )
	print( "read", fileStr2 )
	print( outDict )
	writeOutput( outPre + '.tsv', outDict )
	print( "wrote", outPre + ".tsv" )
	
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print ("python3.4 overlap_regions.py <regions_file_1> <regions_file_2> <out_prefix>")
	else:
		fileStr1 = sys.argv[1]
		fileStr2 = sys.argv[2]
		outPre = sys.argv[3]
		mainFunction( fileStr1, fileStr2, outPre )
		
