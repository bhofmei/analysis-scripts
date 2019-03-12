import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3 gff_to_bed.py [-s=feat-id] <in_gff> <out_bed>

def processInputs(inFileStr, outFileStr, featId):
	if featId.lower() in ['none', 'false']:
		featId = False
	print('Input file: {:s}'.format(inFileStr))
	print('Output file: {:s}'.format(outFileStr))
	readFile( inFileStr, outFileStr, featId )
	print( 'Done' )

def readFile( inFileStr, outFileStr, featId ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in inFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
		# (6) strand (7) frame (8) notes
		chrm = lineAr[0]
		start = int( lineAr[3] ) # this is 1-based
		end = int(lineAr[4] ) # this is 1-based, inclusive
		strand = lineAr[6]
		score = lineAr[5]
		name = '.'
		if featId != False:
			tmp = searchNotes( lineAr[8], featId+'=' )
			if tmp != -1 and tmp != '':
				name = tmp
		# write to output
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		outStr = '{:s}\t{:d}\t{:d}\t{:s}\t{:s}\t{:s}\n'.format( chrm, start-1, end, name, score, strand )
		outFile.write(outStr)
	# end for line
	
	inFile.close()
	outFile.close()
	
def searchNotes( notesStr, searchStr ):
	index = notesStr.find( searchStr )
	if index == -1:
		return ''
	adIndex = index + len( searchStr )
	endIndex = notesStr[adIndex:].find( ';' )
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		newIndex = endIndex+adIndex
		return notesStr[adIndex:newIndex]

def parseInputs( argv ):
	featId = 'ID'
	startInd = 0
	
	for i in range(min(len(argv)-2, 1)):
		if argv[i].startswith('-s='):
			featId = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-'):
			print( 'Options {:s} is not recognzied'.format(argv[i]))
	# end for
	inFileStr = argv[startInd]
	outFileStr = argv[startInd+1]
	
	processInputs(inFileStr, outFileStr, featId)

if __name__ == "__main__":
	if len(sys.argv) != 3 :
		print ("Usage: python gff_to_bed.py [-s=feat-id] <in_gff> <out_bed>")
	else:
		parseInputs( sys.argv[1:] )
