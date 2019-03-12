import sys, math, glob, multiprocessing, subprocess, os
'''
	produces a text file with p-values comparing number of methylated
	positions and unmethylated positions of the genome
	pairwise comparison for all samples include in input file
	input file has header: meth.types line.name ect.
'''

# Usage: pos_methyl_stat.py <in_file> [out_file]

def mainFunction( inFileStr, outPrefix ):
	# read in file
	print( 'Reading in file...' )
	valDict, sampleNamesAr = readInFile( inFileStr )
	# compute pvalues
	print( 'Computing p-values...' )
	pvalDict, outvalDict = computePValues( valDict )
	#print( pvalDict )
	# write output
	outFileStr = outPrefix + '_pval.tsv'
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, outvalDict, sampleNamesAr )
	# write distance matrix type output
	print( 'Writing matrix-like output...' )
	writeMatrixFile( outPrefix, pvalDict, sampleNamesAr  )

def readInFile( inFileStr ):
	outDict = {}
	sampleNamesAr = []
	curMType = None
	inFile = open( inFileStr, 'r' )
	for line in inFile:
		# header
		if line.startswith( 'meth' ) or line.startswith('#'):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) meth.typ (1) line.name (2) num.pos.meth (3) num.pos.unmeth 
		# (4) num.pos.na (5) weighted.meth (6) weighted.unmeth 
		# (7) weighted.overall
		if lineAr[0] != curMType:
			outDict[ lineAr[0] ] = []
			curMType = lineAr[0]
		# regular line
		if lineAr[1] not in sampleNamesAr:
			sampleNamesAr += [ lineAr[1].replace('-','.') ]
		outDict[ lineAr[0] ] += [ ( int(lineAr[2]), int(lineAr[3]) ) ]
	inFile.close()
	return outDict, sampleNamesAr

def computePValues( valDict ):
	pvalDict = {}
	outvalDict = {}
	
	for mType in valDict.keys():
		numLines = len(valDict[mType]) 
		pvalDict[mType] = [ [-1]*numLines for x in range(numLines) ]
		outvalDict[mType] = [ [-1]*numLines for x in range(numLines) ]
		
		# iterate through pairs
		for i in range(numLines):
			for j in range( i ):
				# compute pvalue
				tup1 = valDict[mType][i]
				tup2 = valDict[mType][j]
				p = calculatePVal( tup1, tup2 )
				outvalDict[mType][i][j] = ( tup1[0],tup1[1],tup2[0],tup2[1],p )
				pvalDict[mType][i][j] = p
			# end for j
		# end for i
	# end for mType
	return pvalDict, outvalDict

def calculatePVal( tup1, tup2 ):
	'''
		tup1 and tup2 are tuples of (num meth, num unmethylated) for a sample
		this function calculates the p-value as a two-proprortion z-test
		between these two sets
		return p-value (float)
	'''
	#print( tup1, tup2)
	# P = [(X0) + (X1)]/[(N0)+(N1)]
	# X is reads methylated, N is total reads
	P = float( tup1[0] + tup2[0] ) / float( tup1[1] + tup2[1] )
	p1 = float(tup1[0] ) / tup1[1]
	p2 = float(tup2[0] ) / tup2[1]
	#print( P, p1, p2 )
	# z = (p0) - (p1) / sqrt ( (P*(1-P)/(N0))+ (P*(1-P)/(N1)))
	if P == 0 or P == 1:
		return 1
	z = ( p1 - p2 ) / math.sqrt( (P*(1-P)/float(tup1[1])) + (P*(1-P)/float(tup2[1])) )
	
	# p(x,u,s) = 1/( s *sqrt(2*pi) ) e^-((x-u)^2/(2*s^2))
	return stdNormDist( math.fabs(z) )
	
def stdNormDist( x ):
	'''
		calculate the p-value of two-proportion z test
		calculates 'cumulative probability density' then adjusts
		for 2-tailed and get value of occurring by chance alone
		
		F(x) = 1/2 ( 1 + erf( x/sqrt(2) ) )
		pvalue = 2 * ( 1 - F(x) )
	'''
	cdf = ( 1 + math.erf( x / math.sqrt(2.0) ) ) / 2
	return 2 * ( 1 - cdf )

def writeOutput( outFileStr, pvalDict, sampleNamesAr ):
	
	outFile = open( outFileStr, 'w' )
	# header
	outFile.write( 'mType\tsample1\tsample2\tvalues1\tvalues2\tpvalue\t-log(pval)\n' )
	for mType in pvalDict.keys():
		for i in range(len(sampleNamesAr)):
			for j in range( i ):
				tup = pvalDict[mType][i][j]
				vals1 = '{:d},{:d}'.format(tup[0],tup[1])
				vals2 = '{:d},{:d}'.format(tup[2],tup[3])
				if tup[4] == -1:
					lpval = float( 'nan' )
				else:
					try:
						lpval = -1 * math.log10(tup[4])
					except ValueError:
						lpval = -1
				outStr = '{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:.4f}\t{:.4f}\n'.format( mType, sampleNamesAr[i],sampleNamesAr[j],vals1,vals2,tup[4],lpval )
				outFile.write( outStr )
	# end for mType
	outFile.close()
	
def writeMatrixFile( baseName, pvalDict, sampleNamesAr ):
	
	for mType in pvalDict.keys():
		outFile = open( baseName+'_'+mType+'.tsv', 'w' )
		header = '\t' + '\t'.join( sampleNamesAr )
		outFile.write( header + '\n' )
		for i in range(len(sampleNamesAr)):
			pvals = pvalDict[mType][i]
			outStr = sampleNamesAr[i]
			for k in range(len(pvals)):
				if pvals[k] == -1:
					outStr += '\tNA'
				else:
					try:
						lval = -1 * math.log10(pvals[k])
					except ValueError:
						lval = -1
					outStr += '\t{:.4f}'.format( lval )
			outFile.write( outStr + '\n' )
		outFile.close()
	# end for mType
		

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: pos_methyl_stat.py <in_file> [out_prefix]")
	else:
		inFileStr = sys.argv[1]
		if len( sys.argv ) == 3:
			outPrefix = sys.argv[2]
		else:
			rInd = inFileStr.rfind('.')
			outPrefix= inFileStr[:rInd] + '_out'
		mainFunction( inFileStr, outPrefix )
