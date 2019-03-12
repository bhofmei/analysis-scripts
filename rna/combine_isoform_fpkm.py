import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 combine_isoform_fpkm.py [-o=out_ID] <in_file>

def mainFunction( inFileStr, outID ):
	
	if outID == None:
		outID = 'combined'
	rInd = inFileStr.rfind( '.' )
	outFileStr = inFileStr[:rInd] + '_' + outID + inFileStr[rInd:]
	
	print( 'Processing {:s}...'.format( inFileStr ) )
	readFile( inFileStr, outFileStr )
	print( 'Output written to {:s}'.format( outFileStr ) )

def readFile( inFileStr, outFileStr ):
	
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	isHeader = True
	addIndexAr = []
	
	curGene = None
	tmpLineAr = None
	
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		if isHeader:
			isHeader = False
			for j in range(len(lineAr ) ):
				if lineAr[j].endswith( '_FPKM' ) or lineAr[j].endswith( '_conf_lo' ) or lineAr[j].endswith( '_conf_hi' ):
					addIndexAr += [ j ]
			outFile.write( line )
		else:
			# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id
			# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage
			# (9) S1_FPKM (10) S1_conf_lo (11) S1_conf_hi (12) S1_status (13)+
			# handle non-genes
			if lineAr[1] == '-' or lineAr[1] == 'x':
				curGene = None
				tmpLineAr = None
				outFile.write(line)
			# handle gene when no current gene
			elif curGene == None:
				curGene = lineAr[4]
				tmpLineAr = lineAr
			# handle gene -> genes don't match
			elif curGene != lineAr[4]:
				outFile.write( '\t'.join( tmpLineAr ) + '\n' )
				curGene = lineAr[4]
				tmpLineAr = lineAr
			# handle gene -> genes match
			else:
				tmpLineAr = combineLines( lineAr, tmpLineAr, addIndexAr )
	# end for
	if tmpLineAr != None:
		outFile.write( '\t'.join( tmpLineAr ) + '\n' )
	outFile.close()
	inFile.close()

def combineLines( newLineAr, outLineAr, addIndexAr ):
	# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id
	# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage
	# (9) S1_FPKM (10) S1_conf_lo (11) S1_conf_hi (12) S1_status (13)+
	
	outAr = outLineAr
	# comma-check columns: 0,2,3,5
	for i in [ 0, 2, 3, 5]:
		if newLineAr[i] not in outLineAr[i]:
			outAr[i] += ',' + newLineAr[i]
	# replace with '-' lines: 7
	outAr[7] = '-'
	# add value columns
	for i in addIndexAr:
		outAr[i] =  str(float( newLineAr[i] ) + float( outLineAr[i] ))
	return outAr

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3.4 combine_isoform_fpkm.py [-o=out_ID] <in_file>")
	else:
		outID = None
		sInd = 1
		if sys.argv[1].startswith( '-o=' ):
			outID = sys.argv[1][3:]
			sInd += 1
		inFileStr = sys.argv[sInd]
		mainFunction( inFileStr, outID )
