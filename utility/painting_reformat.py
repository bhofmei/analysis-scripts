import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 painting_reformat.py <in_file> <out_file>

def readFile( inFileStr, outFileStr ):
	
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	header = 'sample\tposition\tvalue\n'
	outFile.write( header )
	prevEnd = 0
	prevSample = None
	for line in inFile:
		if line.startswith( "[" ):
			continue
		lineAr = line.rstrip().split( ' ' )
		samp = lineAr[0]
		if samp != prevSample:
			prevEnd = 0
			prevSample = samp
		for j in range( 1, len(lineAr), 3 ):
			#val = int(lineAr[j])
			val = float(lineAr[j])
			pos1 = int(lineAr[j+1])
			pos2 = int(lineAr[j+2])
			if prevEnd + 1 != pos1:
				outFile.write('{:s}\t{:d}\tNA\n{:s}\t{:d}\tNA\n'.format( samp, prevEnd + 1, samp, pos1 - 1 ) )
			outFile.write( '{:s}\t{:d}\t{:f}\n{:s}\t{:d}\t{:f}\n'.format( samp, pos1, val, samp, pos2, val ) )
			prevEnd = pos2
		# for j
	# for line
	
	outFile.close()
	inFile.close()

if __name__ == "__main__":
	if len(sys.argv) != 3 :
		print ("Usage: python3.4 painting_reformat.py <in_file> <out_file>")
	else:
		inFileStr = sys.argv[1]
		outFileStr = sys.argv[2]
		readFile( inFileStr, outFileStr )
