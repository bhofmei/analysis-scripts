import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 concat_fpkm.py <out_prefix> <fpkm_file1> [fpkm_fileN]*

def processInputs( outPre, fpkmFiles ):
	
	outFileStr = outPre + '_genes.fpkm_tracking'
	fpkmDict = {}
	
	for i in range( len(fpkmFiles) ):
		print( 'Reading {:s}...'.format( fpkmFiles[i] ) )
		indivDict = readFPKMFile( fpkmFiles[i] )
		fpkmDict = addToFPKMDict( fpkmDict, indivDict, i )
	
	sampleNames = getSampleNames( fpkmFiles )
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, fpkmDict, sampleNames )
	print( 'Done.' )


def getSampleNames( fpkmFileAr ):
	sampleAr = []
	for fpkmFile in fpkmFileAr:
		rind = fpkmFile.rfind( '/' )
		n = fpkmFile[rind+1:]
		sampleAr += [ n.replace('genes_','').replace('.fpkm_tracking','') ]
	return sampleAr

def readFPKMFile( fpkmFileStr ):
	fpkmDict = {}
	fpkmFile = open( fpkmFileStr, 'r' )
	
	for line in fpkmFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) tracking_id (1) class_code (2) nearest_ref (3) gene_id 
		# (4) gene_short_id (5) tss_id (6) locus (7) length (8) coverage 
		# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
		if lineAr[12] in ['OK','HIDATA','LOWDATA'] and lineAr[4] != '-':
			fpkmDict[ lineAr[4] ] = ( float(lineAr[9]), lineAr[12], lineAr[6] )
	fpkmFile.close()
	return fpkmDict

def addToFPKMDict( fpkmDict, indivDict, sampleNum ):
	
	# iterate through fpkmDict
	for key in fpkmDict.keys():
		tup = indivDict.get( key )
		
		# not there
		if tup == None:
			fpkmDict[key] += [ (-1, 'NONE') ]
		else:
			fpkmDict[key] += [ (tup[0],tup[1]) ]
			del indivDict[key]
	
	# iterate through rest of genes in indivDict
	for key in indivDict.keys():
		tup = indivDict[key]
		newAr = [ tup[2] ]
		newAr += [ (-1, 'NONE') ] * sampleNum
		newAr += [ ( tup[0], tup[1] ) ]
		fpkmDict[key] = newAr
		
	return fpkmDict

def writeOutput( outFileStr, fpkmDict, sampleNames ):
	
	outFile = open( outFileStr, 'w' )
	header = 'gene_short_name\tlocus'
	for s in sampleNames:
		header += '\t{0:s} FPKM\t{0:s} STATUS'.format( s )
	outFile.write( header + '\n' )
	
	# iterate through genes
	for gene in sorted( fpkmDict.keys() ):
		geneAr = fpkmDict[ gene ]
		outStr = '{:s}\t{:s}'.format( gene, geneAr[0] )
		for i in range(1, len(geneAr) ):
			outStr += '\t{:.4}\t{:s}'.format( geneAr[i][0], geneAr[i][1] )
		outFile.write( outStr + '\n' )
	outFile.close()

def parseInputs( argv ):
	
	outPre = argv[0]
	fpkmFiles = []
	for i in range( 1, len(argv) ):
		fpkmFiles += [ argv[i] ]
	
	processInputs( outPre, fpkmFiles )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 concat_fpkm.py <out_prefix> <fpkm_file1> [fpkm_fileN]*")
	else:
		parseInputs( sys.argv[1:] )
