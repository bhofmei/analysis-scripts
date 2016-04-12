import sys

# Usage: python3.4 gff_description_fpkm_tracking.py <gff_file> <fpkm_file> [fpkm_file]*

def readGFF( gffFileStr ):
	'''
		returns dictionary of genes with their descriptions
		key is gene name and value is description
	'''
	gffFile = open( gffFileStr, 'r' )
	gffDict = {}
	
	for line in gffFile:
		line = line.rstrip()
		# header
		if line.startswith( '#' ):
			continue
		lineAr = line.split( '\t' )
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		
		if lineAr[2] == "gene":
			name = getGeneName( lineAr[8] )
			if name == None:
				continue
			des = getGeneDescription( lineAr[8] )
			gffDict[name] = des
	gffFile.close()
	print( 'GFF file read' )
	return gffDict

def getGeneName( notesStr ):
	search = 'Name='
	index = notesStr.find( search )
	if index == -1:
		return None
	adIndex = index + len( search )
	endIndex = notesStr[adIndex:].find( ';' )
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return notesStr[adIndex:endIndex+adIndex]

def getGeneDescription( notesStr ):
	search = 'description='
	index = notesStr.find( search )
	if index == -1:
		return 'NA'
	adIndex = index + len( search )
	endIndex = notesStr[adIndex:].find( ';' )
	if endIndex == -1:
		return cleanDescription( notesStr[adIndex:] )
	else:
		return cleanDescription( notesStr[adIndex:endIndex+adIndex] )

def cleanDescription( descStr ):
	'''
		variant+surface+glycoprotein+%28VSG%29+%28VSG+427-9%29 ->
		variant surface glycoproten (VSG) (VSG 427-9)
	'''
	# replace '+' with spaces
	descStr1 = descStr.replace( '+', ' ' )
	
	while '%' in descStr1:
		ind = descStr1.find( '%' )
		hxStr = descStr1[ind+1:ind+3]
		hxNum = int( hxStr, 16 )
		rp = chr( hxNum )
		#rp1 = '%' + hxStr
		descStr1 = descStr1.replace( '%' + hxStr, rp )
	return descStr1
		
def readFPKM( fpkmFileStr, outFileStr, gffDict ):
	
	fpkmFile = open( fpkmFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	isHeader = True
	insertIndex = -1
	
	for line in fpkmFile:
		line = line.rstrip()
		lineAr = line.split( '\t' )
		# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
		# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
		# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
		
		if isHeader:
			insertIndex = lineAr.index( 'gene_short_name' )
			isHeader = False
			# write header
			outHeader = '\t'.join( lineAr[:insertIndex+1] ) + '\tgene_description\t' + '\t'.join( lineAr[insertIndex+1:] )
			outFile.write( outHeader + '\n' )
		
		else:
			geneName = lineAr[ insertIndex ]
			geneDes = gffDict.get( geneName )
			if geneDes == None:
				geneDes = 'NA'
			outStr = '\t'.join( lineAr[:insertIndex+1] ) + '\t' + geneDes + '\t' + '\t'.join( lineAr[insertIndex+1:] )
			outFile.write( outStr + '\n' )
	fpkmFile.close()
	outFile.close()
	
def getOutFileStr( fpkmFileStr ):
	
	rind = fpkmFileStr.rfind('.')
	return fpkmFileStr[:rind] + '_desc' + fpkmFileStr[rind:]

def mainFunction( gffFileStr, fpkmFileStrAr ):
	
	gffDict = readGFF( gffFileStr )
	
	for fpkmFileStr in fpkmFileStrAr:
		outFileStr = getOutFileStr( fpkmFileStr )
		readFPKM( fpkmFileStr, outFileStr, gffDict )
		print( 'Finished with', fpkmFileStr )
	print( 'Done.' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 gff_description_fpkm_tracking.py <gff_file> <fpkm_file> [fpkm_file]*")
	else:
		gffFileStr = sys.argv[1]
		fpkmFileStrAr = []
		for i in range( 2, len(sys.argv) ):
			fpkmFileStrAr += [ sys.argv[i] ]
		mainFunction( gffFileStr, fpkmFileStrAr )
