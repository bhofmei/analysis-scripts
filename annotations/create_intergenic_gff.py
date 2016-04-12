import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 create_intergenic_gff.py <gff_file> <out_file>

def readGFF( gffFileStr ):
	
	prevPos = 1
	curChrm = None
	outStr = ''
	
	gffFile = open( gffFileStr, 'r' )
	
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		
		# ignore non-genes
		if lineAr[2] != 'gene':
			continue
		
		chrm = lineAr[0]
		adStart = int( lineAr[3] ) - 100
		adEnd = int( lineAr[4] ) + 100
		
		# check new chrm
		if chrm != curChrm:
			#outStr += '{:s}\t{:d}\t{:d}\n'.format( chrm, 1, adStart )
			prevPos = adEnd
			curChrm = chrm
		# same chrm
		else:
			outStr += '{:s}\t{:d}\t{:d}\n'.format( chrm, min(prevPos, adStart), max(prevPos, adStart) )
			prevPos = adEnd
	gffFile.close()
	return outStr

def mainFunction( gffFileStr, outFileStr ):
	
	outStr = readGFF( gffFileStr )
	outFile = open( outFileStr, 'w' )
	outFile.write( outStr )
	outFile.close()
	print( 'Done.' )

if __name__ == "__main__":
	if len(sys.argv) != 3 :
		print ("Usage: python3.4 create_intergenic_gff.py <gff_file> <out_file>")
	else:
		gffFileStr = sys.argv[1]
		outFileStr = sys.argv[2]
		mainFunction( gffFileStr, outFileStr )
