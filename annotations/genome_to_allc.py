import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: genome_to_allc.py [-o=out_id] [-p=num_proc] <fasta_file>
NUMPROC = 1

def processInputs( fastaFileStr, outID, numProc ):
	# read fasta
	print( 'Reading fasta' )
	fastaFile = FileFASTA( fastaFileStr )
	#fastaDict = fastaFile.getFastaDict( isPickle=False )
	fastaDict = fastaFile.readFasta()
	
	if outID == None:
		outID = fastaFile.fbBasename()
	#outFileStr = 'allc_' + outID + '.tsv'
	
	print( 'Processing with {:d} processor'.format( numProc) )
	pool = multiprocessing.Pool( processes = numProc )
	results = [ pool.apply_async(processChrm, args=(fastaDict[x], x, outID) ) for x in fastaDict.keys() ]
	outputStrings = [ p.get() for p in results ]
	

	print( 'Done' )
	#print( 'Writing output to', outFileStr )
	#with open( outFileStr, 'w' ) as f:
	#	for x in outputStrings:
	#		f.write( x )
	# end with

def processChrm( chrmAr, chrmName, outID ):
	print( ' ' + chrmName )
	outFileStr = 'allc_' + outID + '_'+ chrmName +'.tsv'
	outFile = open( outFileStr, 'w' )
	
	# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
	# (6) methylated
	# note: allc files are 1-index
	
	for i in range(1,len(chrmAr)-2):
		# check forward C
		if chrmAr[i] == 'C':
			pos = i
			context = ''.join(chrmAr[i:i+3])
			outAr = [ chrmName, str(pos), '+', context, '1', '1', '1' ]
			outStr = '\t'.join( outAr ) + '\n'
			outFile.write( outStr )
		# check reverse G
		if chrmAr[i+2] == 'G':
			pos = i+2
			context = reverseComp( chrmAr[i:i+3] )
			outAr = [ chrmName, str(pos), '-', context, '1', '1', '1' ]
			outStr = '\t'.join( outAr ) + '\n'
			outFile.write( outStr )
	# end for i
	outFile.close()
	return outStr
			
def reverseComp( inStr ):
	outStr = ''
	n = len(inStr)
	for i in range(n-1,-1,-1):
		outStr += baseComplement( inStr[i] )
	return outStr

def baseComplement( base ):
	baseAr = ['A','C','G','T']
	revAr = ['T','G','C','A']
	try:
		ind = baseAr.index( base.upper() )
		return revAr[ind]
	except ValueError:
		return 'X'
	

def parseInputs( argv ):
	outID = None
	numProc = NUMPROC
	startInd = 0
	
	for i in range(min(2,len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer...using default', NUMPROC )
				numProc = NUMPROC
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
			
	# end for
	fastaFileStr = argv[startInd]
	processInputs( fastaFileStr, outID, numProc )


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print ("Usage: genome_to_allc.py [-o=out_id] [-p=num_proc] <fasta_file>")
	else:
		parseInputs( sys.argv[1:] )
