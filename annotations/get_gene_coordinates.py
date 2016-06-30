import sys, os
from bioFiles import *

# Usage: python3 get_gene_coordinates.py [-a] [-c=col_num] [-o=out_id] <in_file> <gff_file>

def processInputs( inFileStr, gffFileStr, colNum, outID, isAppend ):
	print( 'Input file:', os.path.basename( inFileStr ) )
	print( 'GFF file:', os.path.basename( gffFileStr ) )
	print( 'Reading gff file' )
	gffFile = FileGFF( gffFileStr )
	gffDict = gffFile.getGeneDict(includeTE=True)
	
	if outID != None:
		outFileStr = '{:s}_coord.tsv'.format( outID )
	else:
		b = os.path.basename( inFileStr )
		rInd = b.rfind( '.' )
		outFileStr = b[:rInd] + '_coord.tsv'
	print( 'Writing output to', outFileStr )
	# read input file and write output
	readFile( inFileStr, outFileStr, gffDict, colNum, isAppend )

def readFile( inFileStr, outFileStr, gffDict, colNum, isAppend ):

	isCSV = ( inFileStr.endswith( '.csv' ) )
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	# header for not append
	if not isAppend:
		outFile.write( '#chrm\tstart\tend\tname\n' )
	
	for line in inFile:
		if isCSV:
			lineAr = line.rstrip().split( ',' )
		else:
			lineAr = line.rstrip().split( '\t' )
		if line.startswith( '#' ):
			if isAppend:
				outFile.write( '{:s}\tgChrm\tgStart\tgEnd\n'.format( '\t'.join( lineAr ) ) )
			continue
		# get gene name
		gName = formatGeneName( lineAr[colNum] )
		gCoord = gffDict.get( gName )
		if gCoord == None:
			print( gName, 'not found in GFF' )
			if isAppend:
				outFile.write( '{:s}\tNA\tNA\tNA\n'.format( '\t'.join( lineAr ) ) )
		else:
			# gCoord is chrm, start, end strand
			gChrm, gStart, gEnd, gStrand = gCoord
			if isAppend:
				outFile.write( '{:s}\t{:s}\t{:d}\t{:d}\n'.format( '\t'.join( lineAr ), gChrm, gStart, gEnd ) )
			else:
				outFile.write( '{:s}\t{:d}\t{:d}\t{:s}\n'.format( gChrm, gStart, gEnd, gName ) )
	# end for line
	inFile.close()
	outFile.close()
	
	
def formatGeneName( inName ):
	rInd = inName.rfind( '.' )
	if rInd != -1:
		return inName[:rInd]
	return inName

def parseInputs( argv ):
	colNum = 0
	outID = None
	isAppend = False
	startInd = 0
	
	for i in range(min(3,len(argv)-2)):
		if argv[i] == '-a':
			isAppend = True
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			try:
				colNum = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: column index must be integer' )
				exit()
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	inFileStr = argv[startInd]
	gffFileStr = argv[startInd+1]
	processInputs( inFileStr, gffFileStr, colNum, outID, isAppend )

def printHelp():
	print ("Usage: python3 get_gene_coordinates.py [-a] [-c=col_num] [-o=out_id] <in_file> <gff_file>")
	print( 'Required:' )
	print( 'in_file\t\ttab-delimited file with gene names' )
	print( 'gff_file\tGFF formatted file' )
	print( 'Optional:' )
	print( '-a\t\tappend coordinates to input file' )
	print( '-o=out_id\toutput identifier [default uses in file name]' )
	print( '-c=col_num\tzero-indexed column index for gene name [default 0]' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
