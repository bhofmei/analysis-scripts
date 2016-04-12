import sys, math, glob, multiprocessing, subprocess, os, random
from bisect import *

# Usage: python3 genes_by_ptu.py [-n=num_genes] <ptu_file> <gff_file>
NUMGENES=1

def processInputs( ptuFileStr, gffFileStr, numGenes ):
	
	info= '#from_script: genes_by_ptu.py; ptu_file: {:s}; gff_file: {:s}; num_genes: {:d}'.format( os.path.basename(ptuFileStr), os.path.basename(gffFileStr), numGenes )
	print( 'Reading GFF file' )
	gffDict = readGFF( gffFileStr )
	print( 'Reading PTU file' )
	ptuDict = readPTU( ptuFileStr )
	
	outFileStr = 'out_ptu_genes_{:d}.tsv'.format( numGenes )
	print( 'Analyzing' )
	outDict = processPTUs( gffDict, ptuDict, numGenes )
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outDict, info )
	print( 'Done' )

def readGFF( gffFileStr ):
	gffDict = {}
	gffFile = open( gffFileStr, 'r' )
	
	for line in gffFile:
		if line.startswith('#'):
				continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) source (2) type (3) start (4) end (5) ?
		# (6) strand (7) ? (8) notes 
		chrmStrand = lineAr[0] + '\t'+lineAr[6]
		if gffDict.get( chrmStrand ) == None:
			gffDict[ chrmStrand ] = []
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		if lineAr[2] == 'gene':
			name = getGeneName( lineAr[8] )
			insort( gffDict[chrmStrand], (start, end, name ) )
	# end for line
	gffFile.close()
	return gffDict

def getGeneName( notesStr):
	search = "Name="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return  notesStr[adIndex:endIndex+adIndex]

def readPTU( ptuFileStr ):
	# return dict chrmStrand: (pos, ptuNum)
	ptuFile = open( ptuFileStr, 'r' )
	ptuDict = {}
	
	for line in ptuFile:
		lineAr = line.rstrip().split('\t')
		# (0) ptu num (1) chrm (2) pos (3) strand
		chrmStrand = lineAr[1] + '\t'+lineAr[3]
		if ptuDict.get( chrmStrand ) == None:
			ptuDict[ chrmStrand ] = []
		ptuDict[ chrmStrand ] += [ ( int(lineAr[2]), int(lineAr[0]) ) ]
	# end for line
	ptuFile.close()
	return ptuDict

def processPTUs( gffDict, ptuDict, numGenes ):
	outDict = {}
	
	# loop through stranded chromosomes
	for chrmStrand in ptuDict.keys():
		# get strandInfo
		isPlus = chrmStrand.endswith( '+' )
		ptuAr = ptuDict[ chrmStrand ]
		gffAr = gffDict.get(chrmStrand)
		# genes are tuples (start, end, name)
		if gffAr != None:	# there are genes on chrm strand
			keys = [ p[0] for p in gffAr ]
			for ptu in ptuAr:
				pos, num = ptu
				outAr = []
				if isPlus:	# plus strand
					genes = getLT( keys, pos, numGenes )
				else:	# negative strand
					genes = getGE( keys, pos, numGenes )
				for i in range(len(genes)):
					gn = gffAr[ genes[i] ][2]
					# ptu, chrm, ptu-pos, gene_dist, gene_name
					outAr += [ (str(num), chrmStrand, str(pos), str(i), gn) ]
				# save to outDict
				outDict[ num ] = outAr
		else:	# there are not genes on the strand
			for ptu in ptuAr:
				pos, num = ptu
				outDict[ num ] = [(num, chrmStrand, str(pos), 'NA', 'NA' )]
	# end for chrmStrand
	return outDict
			
def getLT( a, x, n ):
	'''
		a -> in order array, x -> search key, n-> number to return
	'''
	i = bisect_left(a, x)	# rightmost endpoint
	print( i)
	if i:
		j = max( 0, i - n )	# leftmost endpoint or 0
		l = list(range(j,i))
		l.reverse()
		return l
	return []
	
def getGE( a, x, n ):
	'''
		a -> in order array, x -> search key, n-> number to return
	'''
	i = bisect_left(a, x)	#leftmost endpoint
	if i != len(a):
		j = min( i+n, len(a) )
		l = list(range(i,j))
		return l
	return []

def writeOutput( outFileStr, outDict, info ):
	header = 'ptu_num\tchrm\tstrand\tpos\tgene_dist\tgene_name\n'
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\n' + header )
	
	for num in sorted( outDict.keys() ):
		# gene array
		geneAr = outDict[num]
		for tup in geneAr:
			outFile.write( '\t'.join( tup ) + '\n' )
	# end for num
	outFile.close()
	
def parseInputs( argv ):
	numGenes = NUMGENES
	startInd = 0
	
	for i in range(1,len(argv)):
		if argv[i].startswith( '-n='):
			try:
				numGenes = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of genes must be an integer' )
				exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid argument'.format( argv[i] ) )
			exit()
	# end for
	ptuFileStr = argv[startInd]
	gffFileStr = argv[startInd + 1]
	processInputs( ptuFileStr, gffFileStr, numGenes )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3 genes_by_ptu.py [-n=num_genes] <ptu_file> <gff_file>")
	else:
		parseInputs( sys.argv[1:] )
