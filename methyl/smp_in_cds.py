import sys, math, glob, multiprocessing, subprocess, os, bisect

# Usage: python3.4 smp_in_gene_body.py <smp_file> <gff_file> <gBM_file>

def readGBM( gbmFileStr ):
	
	gbmFile = open( gbmFileStr, 'r' )
	gbmAr = []
	
	for line in gbmFile:
		if line.startswith( 't' ):
			continue
		line = line.rstrip()
		rind = line.rfind( '.' )
		geneName = line[:rind] + '.1'
		gbmAr +=[ geneName ]
	gbmFile.close()
	
	return gbmAr

def readSMP( smpFileStr ):
	
	smpFile = open( smpFileStr, 'r' )
	smpDict = {}
	smpCount = 0
	
	for line in smpFile:
		if line.startswith( '#' ) or line == '\n':
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
	
	return smpDict, smpCount
	
def bisectIndex( a, x ):
	i = bisect.bisect_left( a, x )
	if i != len( a ) and a[i] == x:
		return i
	else:
		return None

def readGFF( gffFileStr, gbmAr ):
	
	gffFile = open( gffFileStr, 'r' )
	gffArCDS = []
	
	for line in gffFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		start = int(lineAr[3])
		end = int(lineAr[4])
		scaffold = lineAr[0]
		strand = lineAr[6]
		#print ( "GFF:",lineAr)
		
		if lineAr[2] == "CDS":
			name, isPrimary = getGeneNameCDS( lineAr[8] )
			if isPrimary and name in gbmAr:
				gffArCDS += [ (scaffold, start, end) ]
	
	gffFile.close()
	return gffArCDS

def getGeneNameCDS (notesStr):
	search = "Parent="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(',')
	name = notesStr[adIndex:endIndex+adIndex]
	isPrimary = name.endswith( '.1' )
	return name, isPrimary

def countSMPs( gffAr, smpDict, outFileStr ):
	
	if outFileStr != None:
		outFile = open( outFileStr, 'w' )
	countGBM = 0
	
	# loop through genes
	for gene in gffAr:
		chrm = gene[0]
		start = gene[1]
		end = gene[2]
		try:
			smpAr = smpDict[ chrm ]
		except KeyError:
			continue
		for i in range( start, end+1 ):
			isIn = bisectIndex( smpAr, i )
			if isIn != None:
				countGBM += 1
				if outFileStr != None:
					outFile.write( '{:s}\t{:d}\n'.format(chrm,i) )
	if outFileStr != None:
		outFile.close()
	return countGBM

def processInputs( smpFileStr, gffFileStr, gbmFileStr ):
	outFileStr = smpFileStr[:(smpFileStr.rfind('.'))] + '_gbm_cds.tsv'
	print( 'Reading SMP file...' )
	smpDict, smpCount = readSMP( smpFileStr )
	print( 'Reading gBM file...' )
	gbmAr = readGBM( gbmFileStr )
	print( 'Reading GFF file...' )
	gffArGBM = readGFF( gffFileStr, gbmAr )
	print( 'Sorting SMPs...' )
	countGBM = countSMPs( gffArGBM, smpDict, outFileStr )
	print( '\n{:d} CG-SMPs in gBM genes\n{:d} total CG-SMPs\n'.format( countGBM, smpCount ) )
	
if __name__ == "__main__":
	if len(sys.argv) != 4 :
		print ("Usage: python3.4 smp_in_gene_body.py <smp_file> <gff_file> <gBM_file>")
	else:
		smpFileStr = sys.argv[1]
		gffFileStr = sys.argv[2]
		gbmFileStr = sys.argv[3]
		processInputs( smpFileStr, gffFileStr, gbmFileStr )
