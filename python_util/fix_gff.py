import copy

def fixGrapeGFF( inFileStr, outFileStr ):
	# missing counts
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	count = 0
	for line in inFile:
		# check (2) for feature types
		lineAr = line.rstrip().split('\t')
		if lineAr[2] == 'gene':
			count = 0
		elif lineAr[2] == 'mRNA':
			count += 1
			# (8) has notes
			lineAr[8] = updateID( lineAr[8], 'ID=', ';', count )
		elif lineAr[2] in ['CDS','three_prime_UTR', 'five_prime_UTR']:
			lineAr[8] = updateID( lineAr[8],'ID=', '.', count )
			lineAr[8] = updateID( lineAr[8], 'Parent=', ';', count )
		# write
		outFile.write( '{:s}\n'.format( '\t'.join(lineAr) ) )
	# end for
	inFile.close()
	outFile.close()

def fixBeetGFF( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	#5-utr, cds, 3-utr
	count = [1]*3
	for line in inFile:
		# check for feature types
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		if lineAr[2] not in ['gene', 'mRNA', 'five_prime_UTR', 'CDS', 'three_prime_UTR' ]:
			continue
		elif lineAr[2] == 'gene':
			count = [1]*3
		elif lineAr[2] == 'mRNA':
			count = [1]*3
		# types to count
		elif lineAr[2] == 'five_prime_UTR':
			lineAr[8] = updateUTR( lineAr[8], count[0], 'five_prime_UTR' )
			count[0] += 1
		elif lineAr[2] == 'CDS':
			lineAr[8] = updateID( lineAr[8], 'ID=', ';', count[1] )
			count[1] += 1
		elif lineAr[2] == 'three_prime_UTR':
			lineAr[8] = updateUTR( lineAr[8], count[2], 'three_prime_UTR' )
			count[2] += 1
		outFile.write( '{:s}\n'.format( '\t'.join(lineAr) ) )
	# end for line
	inFile.close()
	outFile.close()

def fixYeast( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	undefCount = 1
	for line in inFile:
		# all lines should be genes but we'll check do some checks anyways
		lineAr = line.rstrip().split( '\t' )
		if line.startswith( '#' ) or len( lineAr ) < 9 or lineAr[2] != 'gene':
			continue
		elif lineAr[8] == 'UNDEF':
			lineAr[8], id, name = updateUndef( undefCount )
			undefCount += 1
		elif ',' not in lineAr[8]:
			lineAr[8], id, name = updatePoorAnnotated( lineAr[8] )
		else:
			lineAr[8], id, name = updateYeastGene( lineAr[8] )
		newLines = createFeatures( lineAr[:8], id, name )
		outFile.write( '{:s}\n{:s}'.format( '\t'.join( lineAr ), newLines ) )
	# end for line
	inFile.close()
	outFile.close()

def fixMaize( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		if line.startswith( '#' ) or len( lineAr ) < 9:
			continue
		# check for "ensembl" label
		if lineAr[1] == '.':
			lineAr[1] = 'ensembl'
		# mRNA-remove mRNA from ID
		if lineAr[2] == 'mRNA':
			lineAr[8] = lineAr[8].replace('mRNA:','')
		# CDS-remove CDS from ID, add CDS # from rank, remove mRNA from parent
		elif lineAr[2] == 'CDS':
			lineAr[8] = updateMaizeCDS( lineAr[8] )
		# UTR-add ID
		elif lineAr[2] == 'five_prime_UTR':
			lineAr[8] = updateMaizeUTR( lineAr[8], '5UTR' )
		elif lineAr[2] == 'three_prime_UTR':
			lineAr[8] = updateMaizeUTR( lineAr[8], '3UTR' )
		# write
		outFile.write( '{:s}\n'.format( '\t'.join(lineAr) ) )
	# end for line
	inFile.close()
	outFile.close()

def updateMaizeCDS( notesStr ):
	# replace "mRNA" and "CDS"
	newNotes = notesStr.replace( 'CDS:', '' ).replace( 'mRNA:', '' )
	# get rank
	index = newNotes.find( 'rank=' )
	if index == -1:
		print( 'ERROR: no rank for', newNotes )
		exit()
	adIndex = index + len( 'rank=' )
	endIndex = newNotes[adIndex:].find( ';' )
	if endIndex == -1:
		rank = newNotes[adIndex:]
	else:
		rank = newNotes[adIndex:adIndex+endIndex]
	iRank = int( rank )
	# add rank to ID
	jndex = newNotes.find( 'ID=' )
	adJndex = jndex + len( 'ID=' )
	endJndex = newNotes[adJndex:].find( ';' )
	if endJndex == -1:
		return newNotes + '_C{:03d}'.format( iRank )
	else:
		newJndex = endJndex+adJndex
		return newNotes[:newJndex] + '_C{:03d}'.format( iRank ) + newNotes[newJndex:]

def updateMaizeUTR( notesStr, utrType ):
	# get ID
	index = notesStr.find( 'Parent=' )
	rindex = notesStr.rfind( '_' )
	geneID = notesStr[ index+len('Parent=') : rindex ]
	utrID = 'ID={:s}_{:s};'.format( geneID, utrType )
	return utrID + notesStr

def updateUndef( count ):
	id = 'UNDEF_{:d}'.format( count )
	outStr = 'ID={:s};Name={:s}'.format( id, id )
	return outStr, id, id

def updatePoorAnnotated( notesStr ):
	id = notesStr
	outStr = 'ID={:s};Name={:s}'.format( id, id )
	return outStr, id, id

def updateYeastGene( notesStr ):
	# split by ; - indicates paralogs
	notesAr = notesStr.split( ';' )
	# get gene id name from first one
	info = notesAr[0].split( ',' )
	# (0) ID (1) SacCer (2) chrm (3) start (4) end (5) name (6) # (7) #
	id = info[0]
	name = info[5]
	notes = '{:s},{:s}'.format( ','.join( info[1:5] ), ','.join( info[6:] ) )
	outStr = 'ID={:s};Name={:s};Notes:{:s}'.format( id, name, notes )
	for i in range(1,len(notesAr)):
		outStr += ';Homolog{:d}={:s}'.format( i, notesAr[i] )
		info = notesAr[i].split( ',' )
		l1,l2 = info[0], info[5]
		outStr += ';Alias={:s}'.format( l1 )
		if l1 != l2:
			outStr += ';Alias={:s}'.format( l2 )
	# end for
	return outStr, id, name

def createFeatures( inAr, id, name ):
	# mRNA
	inAr[2] = 'mRNA'
	notes = 'ID={:s}m;Name={:s}m;Parent={:s}'.format( id, name, id )
	outStr = '{:s}\t{:s}\n'.format( '\t'.join( inAr ), notes )
	inAr[2] = 'CDS'
	notes = notes = 'ID={:s}m.CDS.1;Parent={:s}m'.format( id, id )
	outStr += '{:s}\t{:s}\n'.format( '\t'.join( inAr ), notes )
	return outStr	

def updateUTR( notesStr, count, replacement ):
	searchStr = 'ID='
	notesStr = notesStr.replace( 'UTR', replacement)
	index = notesStr.find( searchStr )
	adIndex = index + len( searchStr )
	endIndex = notesStr[adIndex:].find( ';' )
	if endIndex == -1:
		return notesStr + '.' + str(count)
	else:
		newIndex = endIndex+adIndex
		return notesStr[:newIndex] + '.' + str(count) + notesStr[newIndex:]

def updateID( notesStr, searchStr, endSearchStr, count ):
	index = notesStr.find( searchStr )
	adIndex = index + len( searchStr )
	endIndex = notesStr[adIndex:].find( endSearchStr )
	if endIndex == -1:
		return notesStr + '.' + str(count)
	else:
		newIndex = endIndex+adIndex
		return notesStr[:newIndex] + '.' + str(count) + notesStr[newIndex:]

def getID( notesStr, searchStr ):
	index = notesStr.find( searchStr )
	adIndex = index + len( searchStr )
	endIndex = notesStr[adIndex:].find( ';' )
	if endIndex == -1:
		return notesStr
	else:
		newIndex = endIndex+adIndex
		return notesStr[:newIndex]

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

def fixPeanut( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	for line in inFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		if 'gene' not in lineAr[2]:
			outFile.write( line )
		else:
			note = updatePeanutNote( lineAr[8] )
			lineAr[8] += note
			outFile.write( '\t'.join( lineAr ) + '\n' )
	#end for line
	inFile.close()
	outFile.close()

def updatePeanutNote( notesStr ):
	search = "Description="
	index = notesStr.find(search)
	if index == -1:
		return ""
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		desc = notesStr[adIndex:]
	else:
		desc = notesStr[adIndex:endIndex+adIndex]
	
	# check for %3B
	index = desc.find( '%3B' )
	if index == -1:
		note = desc
	else:
		note = desc[:index]
	return ';Note=' + note 

#### maize v4 ####
	
def fixMaize4( inFileStr, outFileStr, ncOutFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	ncFile = open(ncOutFileStr, 'w')
	
	geneAr = []
	mrnaAr = []
	strand = None
	isNoncoding = False
	for line in inFile:
		line = line.rstrip()
		if line == '###':
			#print( line  )
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				if isNoncoding:
					tmpStr = handleNoncoding_maize( geneAr )
					ncFile.write( tmpStr )
				else:
					tmpStr = handleGene( geneAr, strand )
					outFile.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
			isNoncoding = False
		elif line.startswith( '#' ) :
			continue
		else:
			lineAr = line.split( '\t' )
			if len(lineAr ) < 9:
				continue
			# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
			# (6) strand (7) frame (8) notes 

			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			if tmp['feature'] == 'gene':
				geneAr.append( tmp )
				strand = tmp['strand']
			elif tmp['feature'].endswith( 'gene' ):
				geneAr.append( tmp )
				strand = tmp['strand']
				isNoncoding = True
			elif tmp['feature'] == 'mRNA':
				#print( mrnaAr )
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				mrnaAr = [tmp]
			elif tmp['feature'].endswith( 'RNA' ) or tmp['feature'] == 'transcript':
				if len(mrnaAr) > 0:
					geneAr.append(mrnaAr)
				mrnaAr = [tmp]
			elif (tmp['feature'] in [ 'contig', 'chromosome']) or (tmp['feature'] == 'exon' and not isNoncoding):
				continue
			else:
				mrnaAr.append( tmp )
		# end else
	# end for line
			
	# end for line
	inFile.close()
	outFile.close()
	ncFile.close()
	print( 'Done' )

def handleNoncoding_maize( inAr ):
	gene = inAr[0]
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		mRNA = mrna[0]
		exons = mrna[1:]
		outStr += outputGene( mRNA )
		for exon in exons:
			outStr += outputGene( exon )
	# end for mrna
	return outStr

def handleGene( inAr, strand ):
	#print( inAr )
	gene = inAr[0]
	try:
		gene['notes'] = gene['notes'].replace( 'gene:', '' )
	except TypeError:
		print (gene)
		raise
		exit()
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		if strand == '+':
			outStr += handlePositive( mrna )	
		else:
			outStr += handleNegative( mrna)
	# end for
	return outStr
	
		
def handlePositive( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mRNA['notes'] = mRNA['notes'].replace('gene:', '' ).replace('transcript:', '' )
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		feat['notes'] = feat['notes'].replace('transcript:', '' )
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [feat]
	#print( mrnaID)	
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetPositive( mrnaID, UTR5, CDS, UTR3 )
	
	return outStr
	
def createGeneSetPositive( mrnaID, UTR5, CDS, UTR3 ):
	
	outStr = ''
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = i+1
		#notes = mrnaID+'_five.prime.UTR'+str(j)
		x = UTR5[i]
		#x['notes'] += ';' + notes
		x['notes'] = updateMaize4UTR5( mrnaID, x['notes'], j )
		outStr += outputGene( x )
	
	# loop through CDS
	for i in range(len(CDS)):
		j = i+1
		x = CDS[i]
		x['notes'] = updateMaize4CDS( x['notes'], j )
		outStr += outputGene( x )
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = i+1
		#notes = mrnaID+'_three.prime.UTR'+str(j)
		x = UTR3[i]
		#x['notes'] += ';' + notes
		x['notes'] = updateMaize4UTR3( mrnaID, x['notes'], j )
		outStr += outputGene( x )
	return outStr

def updateMaize4CDS( notes, count ):
	#newNotes = notes.replace('CDS:', '' ).replace(':cds','')
	newNotes = notes.replace(':CDS1', '' )
	jndex = newNotes.find( 'ID=' )
	adJndex = jndex + len( 'ID=' )
	endJndex = newNotes[adJndex:].find( ';' )
	if endJndex == -1:
		return newNotes + ':CDS{:d}'.format( count )
	else:
		newJndex = endJndex+adJndex
		return newNotes[:newJndex] + ':CDS{:d}'.format( count ) + newNotes[newJndex:]
		
def updateMaize4UTR5( mrnaID, notes, count ):
	#newNotes = notes.replace(':five_prime_utr','')
	newNotes = notes.replace(':five_prime_UTR1','')
	jndex = newNotes.find( 'ID=' )
	if jndex == -1:
		#return newNotes + ';{:s}_five.prime.UTR{:d}'.format( mrnaID, count )
		return newNotes + ';{:s}:five_prime_UTR{:d}'.format( mrnaID, count )
	adJndex = jndex + len( 'ID=' )
	endJndex = newNotes[adJndex:].find( ';' )
	if endJndex == -1:
		#return newNotes + '_five.prime.UTR{:d}'.format( count )
		return newNotes + ':five_prime_UTR{:d}'.format( count )
	else:
		newJndex = endJndex+adJndex
		#return newNotes[:newJndex] + '_five.prime.UTR{:d}'.format( count ) + newNotes[newJndex:]
		return newNotes[:newJndex] + ':five_prime_UTR{:d}'.format( count ) + newNotes[newJndex:]
		
def updateMaize4UTR3( mrnaID, notes, count ):
	#newNotes = notes.replace(':three_prime_utr','')
	newNotes = notes.replace(':three_prime_UTR1','')
	jndex = newNotes.find( 'ID=' )
	if jndex == -1:
		#return newNotes + ';{:s}_three.prime.UTR{:d}'.format( mrnaID, count )
		return newNotes + ';{:s}:three_prime_UTR{:d}'.format( mrnaID, count )
	adJndex = jndex + len( 'ID=' )
	endJndex = newNotes[adJndex:].find( ';' )
	if endJndex == -1:
		#return newNotes + '_three.prime.UTR{:d}'.format( count )
		return newNotes + ':three_prime_UTR{:d}'.format( count )
	else:
		newJndex = endJndex+adJndex
		#return newNotes[:newJndex] + '_three.prime.UTR{:d}'.format( count ) + newNotes[newJndex:]
		return newNotes[:newJndex] + ':three_prime_UTR{:d}'.format( count ) + newNotes[newJndex:]
	
def handleNegative( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mRNA['notes'] = mRNA['notes'].replace('gene:', '' ).replace('transcript:', '' )
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		feat['notes'] = feat['notes'].replace('transcript:', '' )
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [ feat ]
		
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetNegative( mrnaID, UTR5, CDS, UTR3 )
	return outStr

def createGeneSetNegative( mrnaID, UTR5, CDS, UTR3 ):
	outStr = ''
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = len(UTR3) - i
		#notes =  mrnaID+'_three.prime.UTR'+str(j) 
		x = UTR3[i]
		#x['notes'] += ';' + notes
		x['notes'] = updateMaize4UTR3( mrnaID, x['notes'], j )
		outStr += outputGene( x )
		
	# loop through CDS
	for i in range(len(CDS)):
		j = len(CDS) - i
		x = CDS[i]
		x['notes'] = updateMaize4CDS( x['notes'], j )
		outStr += outputGene( x )
	
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = len(UTR5) - i
		#notes =  mrnaID+'_five.prime.UTR'+str(j) 
		x = UTR5[i]
		#x['notes'] += ';' + notes
		x['notes'] = updateMaize4UTR5( mrnaID, x['notes'], j )
		outStr += outputGene( x )
		
	#print(outAr)
	return outStr

#### Human ####
def fixHuman( inFileStr, outFileGeneStr, outFileRNAStr ):
	rnaGenes = [ 'miRNA_gene', 'lincRNA_gene', 'rRNA_gene', 'snRNA_gene', 'snoRNA_gene' ]
	rnaTranscripts = ['lincRNA', 'miRNA', 'rRNA', 'snRNA' ,'snoRNA' ]
	geneTypes = ['gene', 'VD_gene_segment', 'C_gene_segment', 'J_gene_segment', 'V_gene_segment', 'processed_transcript', 'mt_gene', 'RNA','pseudogene' ]
	ignoreAr = ['chromosome', 'biological_region', 'supercontig']
	secondaryTranscripts = ['mRNA', 'pseudogene', 'transcript', 'aberrant_processed_transcript', 'VD_gene_segment', 'C_gene_segment', 'J_gene_segment', 'V_gene_segment', 'processed_pseudogene', 'NMD_transcript_variant', 'pseudogenic_transcript','processed_transcript']
	cdsTypes = ['mRNA'] #, 'pseudogene' ]
	
	inFile = open( inFileStr, 'r' )
	outFileGene = open( outFileGeneStr, 'w' )
	outFileRNA = open( outFileRNAStr, 'w' )
	
	geneAr = []
	mrnaAr = []
	strand = None
	type = None
	
	for line in inFile:
		line = line.rstrip()
		if line == '###' :
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				#print( geneAr )
				tmpStr = handleGeneHuman( geneAr, strand )
				if type == 'rna':
					outFileRNA.write( tmpStr )
				else:
					outFileGene.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
			type = None
		elif line.startswith( '#' ):
			continue
		else:
			lineAr = line.split( '\t' )
			if len( lineAr ) < 9:
				continue
			# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
			# (6) strand (7) frame (8) notes
			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			
			feat = tmp['feature']
			# ignore
			if feat in ignoreAr:
				continue
			# gene features
			elif feat in rnaGenes or (feat in geneTypes and type == None):
				strand = tmp['strand']
				geneAr.append( tmp )
				type = 'gene' if feat in geneTypes else 'rna'
			# second level features
			elif feat in secondaryTranscripts or feat in rnaTranscripts:
				if feat in cdsTypes:
					type = 'geneA'
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				#print( mrnaAr )
				mrnaAr = [tmp]
			# exon/CDS
			else:
				if type == 'geneA' and feat == 'exon':
					continue
				mrnaAr.append( tmp )
		# end else
	# end for line
			
def handleGeneHuman( inAr, strand ):
	#print( inAr )
	gene = inAr[0]
	gene['notes'] = gene['notes'].replace( 'gene:', '' )
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		if strand == '+':
			outStr += handlePositiveHuman( mrna )	
		else:
			outStr += handleNegativeHuman( mrna)
	# end for
	return outStr
	
		
def handlePositiveHuman( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mRNA['notes'] = mRNA['notes'].replace('gene:', '' ).replace('transcript:', '' )
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		feat['notes'] = feat['notes'].replace('transcript:', '' )
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS' or feat['feature'] == 'exon':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [feat]
	#print( mrnaID)	
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetPositiveHuman( mrnaID, UTR5, CDS, UTR3 )
	
	return outStr
	
def handleNegativeHuman( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mRNA['notes'] = mRNA['notes'].replace('gene:', '' ).replace('transcript:', '' )
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		feat['notes'] = feat['notes'].replace('transcript:', '' )
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS' or feat['feature'] == 'exon':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [ feat ]
		
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetNegativeHuman( mrnaID, UTR5, CDS, UTR3 )
	return outStr

def updateHumanCDS( notes, mrnaID, count ):
	newNotes = notes.replace('CDS:', '' )
	jndex = newNotes.find( 'ID=' )
	# if it doesn't have an id, it is an exon
	if jndex == -1:
		return 'ID={:s}_EXON{:d};{:s}'.format( mrnaID, count, newNotes )
	adJndex = jndex + len( 'ID=' )
	endJndex = newNotes[adJndex:].find( ';' )
	if endJndex == -1:
		return newNotes + '_CDS{:d}'.format( count )
	else:
		newJndex = endJndex+adJndex
		return newNotes[:newJndex] + '_CDS{:d}'.format( count ) + newNotes[newJndex:]

def createGeneSetPositiveHuman( mrnaID, UTR5, CDS, UTR3 ):
	
	outStr = ''
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = i+1
		notes = mrnaID+'_five.prime.UTR'+str(j)
		x = UTR5[i]
		x['notes'] += ';' + notes
		outStr += outputGene( x )
	
	# loop through CDS
	for i in range(len(CDS)):
		j = i+1
		x = CDS[i]
		x['notes'] = updateHumanCDS( x['notes'], mrnaID, j )
		outStr += outputGene( x )
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = i+1
		notes = mrnaID+'_three.prime.UTR'+str(j)
		x = UTR3[i]
		x['notes'] += ';' + notes
		outStr += outputGene( x )
	return outStr

def createGeneSetNegativeHuman( mrnaID, UTR5, CDS, UTR3 ):
	outStr = ''
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = len(UTR3) - i
		notes =  mrnaID+'_three.prime.UTR'+str(j) 
		x = UTR3[i]
		x['notes'] += ';' + notes
		outStr += outputGene( x )
		
	# loop through CDS
	for i in range(len(CDS)):
		j = len(CDS) - i
		x = CDS[i]
		x['notes'] = updateHumanCDS( x['notes'], mrnaID, j )
		outStr += outputGene( x )
	
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = len(UTR5) - i
		notes =  mrnaID+'_five.prime.UTR'+str(j) 
		x = UTR5[i]
		x['notes'] += ';' + notes
		outStr += outputGene( x )
		
	#print(outAr)
	return outStr

#### beijerinckii ####
def fix_beijerinckii( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	mrnaAr = []
	for line in inFile:
		line = line.rstrip()
		if line == '###':
			#print( line  )
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				tmpStr = handleGene_beijerinckii( geneAr )
				outFile.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
		elif line.startswith( '#' ) :
			continue
		else:
			lineAr = line.split( '\t' )
			if len(lineAr ) < 9:
				continue
			# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
			# (6) strand (7) frame (8) notes 

			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			if tmp['feature'] == 'gene':
				geneAr.append( tmp )
			elif tmp['feature'] == 'mRNA':
				#print( mrnaAr )
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				mrnaAr = [tmp]
			else:
				mrnaAr.append( tmp )
		# end else
	# end for line
			
	# end for line
	inFile.close()
	outFile.close()
	print( 'Done' )

def handleGene_beijerinckii( inAr ):
	
	gene = inAr[0]
	#print( gene )
	#gene['notes'] = gene['notes'].replace( 'gene:', '' )
	
	mrnas = inAr[1:]
	## no mRNA and no CDS
	if len(mrnas) == 0:
		tmpStr = createRNACDS_beijerinckii( gene )
	# no CDS
	elif len(mrnas) == 1:
		print( '-', mrnas[0][0] )
		tmpStr = createRNA_beijerinckii( gene, mrnas[0][0] )
	else:
		print( 'ERROR:', gene )
		exit()	
	
	outStr = '####\n' + tmpStr
	return outStr

def createRNACDS_beijerinckii( geneInfo ):
	
	outStr = outputGene( geneInfo )
	## mRNA
	mRNAInfo = geneInfo.copy()
	
	geneID = searchNotes( geneInfo['notes'], 'ID=')
	mRNAID = geneID + '_T1'
	cdsID = geneID + '_T1_CDS1'
	geneName = searchNotes( geneInfo['notes'], 'Name=')
	mRNAName = geneName + '_T1'
	cdsName = geneName + '_T1_CDS1'
	print(geneID, geneName)
	isTRNA = ('tRNA' in geneInfo['notes'] )
	mRNAInfo['feature'] = ('tRNA' if isTRNA else 'mRNA')
	mRNAInfo['frame'] = '.'
	mRNAInfo['score'] = '.'
	mRNAInfo['source'] = '.'
	mRNAInfo['notes'] = 'ID={:s};Name={:s};Parent={:s}'.format( mRNAID, mRNAName, geneID )
	#print( '--', mRNAInfo )
	outStr += outputGene( mRNAInfo )
	## CDS
	cdsInfo = geneInfo.copy()
	
	cdsInfo['feature'] = 'CDS'
	cdsInfo['source'] = '.'
	cdsInfo['notes'] = 'ID={:s};Name={:s};Parent={:s}'.format( cdsID, cdsName, mRNAID )
	#print( '---', cdsInfo )
	outStr += outputGene( cdsInfo )
	return outStr
	
def createRNA_beijerinckii( geneInfo, cdsInfo ):
	outStr = outputGene( geneInfo )
	mrnaInfo = geneInfo.copy()
	#print( '--', mrnaInfo )
	geneID = searchNotes( geneInfo['notes'], 'ID=')
	mRNAID = geneID + '_T1'
	geneName = searchNotes( geneInfo['notes'], 'Name=')
	mrnaName = geneName + '_T1'
	mrnaInfo['feature'] = 'mRNA'
	mrnaInfo['frame'] = '.'
	mrnaInfo['score'] = '.'
	mrnaInfo['source'] = '.'
	mrnaInfo['notes'] = 'ID={:s};Name={:s};Parent={:s}'.format( mRNAID, mrnaName, geneID )
	outStr += outputGene( mrnaInfo )
	# fix up CDS
	cdsInfo['notes'] = update_beijerinckii_cds( cdsInfo['notes'], mRNAID )
	outStr += outputGene( cdsInfo )
	return outStr
	
def update_beijerinckii_cds( notes, mrnaID ):
	# original parent
	oP = searchNotes( notes, 'Parent=' )
	return notes.replace(oP, mrnaID )

#### Conringa ####
def fixConringa( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	mrnaAr = []
	strand = None
	for line in inFile:
		line = line.rstrip()
		if line == '###':
			#print( line  )
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				tmpStr = handleGene_conringa( geneAr, strand )
				outFile.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
		elif line.startswith( '#' ) :
			continue
		else:
			lineAr = line.split( '\t' )
			if len(lineAr ) < 9:
				continue
			# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
			# (6) strand (7) frame (8) notes 

			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			if tmp['feature'] == 'gene':
				geneAr.append( tmp )
				strand = tmp['strand']
			elif tmp['feature'] == 'mRNA':
				#print( mrnaAr )
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				mrnaAr = [tmp]
			else:
				mrnaAr.append( tmp )
		# end else
	# end for line
	inFile.close()
	outFile.close()
	print( 'Done' )

def handleGene_conringa( inAr, strand ):
	gene = inAr[0]
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		if strand == '+':
			outStr += handlePositive_conringa( mrna )	
		else:
			outStr += handleNegative_conringa( mrna )
	# end for
	return outStr

def handlePositive_conringa( inAr ):
	mRNA = inAr[0]
	outStr = outputGene( mRNA )
	mrnaID = searchNotes( mRNA['notes'], 'ID=' )
	UTR5=[]
	CDS=[]
	UTR3=[]
	
	cnt = 0
	fStart = -1
	fEnd = -1
	start = -1
	end = -1
	#print( inAr[1:] )
	for feat in inAr[1:]:
		if feat['feature'] == 'start_codon':
			start = feat['start'] - 1
			cnt += 1
			if len(UTR5) > 0:
				tmp = UTR5[-1]
				fEnd = tmp['end']
				tmp['end'] = start
		elif feat['feature'] == 'stop_codon':
			end = feat['end'] + 1
			cnt += 1
			tmp = CDS[-1]
			if tmp['feature'] != 'exon':
				tmp = copy.copy(UTR5[-1])
				tmp['end'] = fEnd
			else:
				CDS = CDS[:-1]
			tmp['start'] = end
			UTR3 = [tmp]	
		elif feat['feature'] == 'exon' and cnt == 0:
			UTR5 += [ feat ]
		elif feat['feature'] == 'exon' and cnt == 1:
			CDS += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'exon' and cnt == 2:
			UTR3 += [ feat ]
	# end for feat
	
	# clean up CDS
	newCDS = []
	for x in CDS:
		if x['feature'] == 'CDS':
			newCDS += [x]
	outStr += createGeneSetPositive_conringa( mrnaID, UTR5, newCDS, UTR3 )
	return outStr

def createGeneSetPositive_conringa( mrnaID, UTR5, CDS, UTR3 ):
	outStr = ''
	
	# loop through UTR5
	for i in range(len(UTR5)):
		j = i+1
		x = UTR5[i]
		x['notes'] = 'ID={:s}.five_prime_UTR.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		x['feature'] = 'five_prime_UTR'
		outStr += outputGene( x )
	# loop through CDS
	for i in range(len(CDS)):
		j = i+1
		x = CDS[i]
		x['notes'] = 'ID={:s}.CDS.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	# loop through UTR3
	for i in range(len(UTR3)):
		j = i+1
		x = UTR3[i]
		x['notes'] = 'ID={:s}.three_prime_UTR.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		x['feature'] = 'three_prime_UTR'
		outStr += outputGene( x )
	return outStr

def handleNegative_conringa( inAr ):
	mRNA = inAr[0]
	outStr = outputGene( mRNA )
	mrnaID = searchNotes( mRNA['notes'], 'ID=' )
	UTR5=[]
	CDS=[]
	UTR3=[]	
	
	cnt = 0
	fStart = -1
	fEnd = -1
	start = -1
	end = -1
	for feat in inAr[1:]:
		#print( cnt, '-', UTR3, '-', CDS, '-', UTR5 )
		#print( feat )
		if feat['feature'] == 'stop_codon':
			start = feat['start'] - 1
			cnt += 1
			if len(UTR3) > 0:
				tmp = UTR3[-1]
				fEnd = tmp['end']
				tmp['end'] = start
		elif feat['feature'] == 'start_codon':
			end = feat['end'] + 1
			cnt += 1
			tmp = CDS[-1]
			if tmp['feature'] != 'exon':
				tmp = copy.copy( UTR3[-1] )
				tmp['end'] = fEnd
			else:
				CDS = CDS[:-1]
			tmp['start'] = end
			UTR5 = [tmp]
			
		elif feat['feature'] == 'exon' and cnt == 0:
			UTR3 += [ feat ]
		elif feat['feature'] == 'exon' and cnt == 1:
			CDS += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'exon' and cnt == 2:
			UTR5 += [ feat ]
	# end for feat
	
	# clean up CDS
	newCDS = []
	for x in CDS:
		if x['feature'] == 'CDS':
			newCDS += [x]
	outStr += createGeneSetNegative_conringa( mrnaID, UTR5, newCDS, UTR3 )
	return outStr

def createGeneSetNegative_conringa( mrnaID, UTR5, CDS, UTR3 ):
	outStr = ''
	
	# loop through UTR3
	for i in range(len(UTR3)):
		j = len(UTR3) - i
		x = UTR3[i]
		x['notes'] = 'ID={:s}.three_prime_UTR.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		x['feature'] = 'three_prime_UTR'
		outStr += outputGene( x )
	# loop through CDS
	for i in range(len(CDS)):
		j = len(CDS) - i
		x = CDS[i]
		x['notes'] = 'ID={:s}.CDS.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	# loop through UTR5
	for i in range(len(UTR5)):
		j = len(UTR5) - i
		x = UTR5[i]
		x['notes'] = 'ID={:s}.five_prime_UTR.{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		x['feature'] = 'five_prime_UTR'
		outStr += outputGene( x )
	return outStr

#### C. elegans ####
def fixCelegans( inFileStr, outFileStr, ncOutFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	ncFile = open( ncOutFileStr, 'w' )
	
	geneAr = []
	mrnaAr = []
	strand = None
	isNoncoding = False
	
	for line in inFile:
		line = line.rstrip()
		if line == '###':
			#print( line  )
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				if isNoncoding:
					tmpStr = handleNoncoding_celegans( geneAr )
					ncFile.write( tmpStr )
				else:
					tmpStr = handleGene_celegans( geneAr, strand )
					outFile.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
			isNoncoding = False
		elif line.startswith( '#' ) :
			continue
		else:
			lineAr = line.split( '\t' )
			if len(lineAr ) < 9:
				continue
			# (0) chrm (1) source (2) feature (3) start (4) stop (5) score
			# (6) strand (7) frame (8) notes 

			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			if tmp['feature'] == 'gene':
				geneAr.append( tmp )
				strand = tmp['strand']
			elif tmp['feature'] in ['mRNA', 'pseudogenic_transcript', 'nc_primary_transcript']:
				#print( mrnaAr )
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				mrnaAr = [tmp]
			elif tmp['feature'].endswith('RNA') or tmp['feature'] == 'miRNA_primary_transcript':
				if len(mrnaAr) > 0:
					geneAr.append(mrnaAr)
				mrnaAr = [tmp]
				isNoncoding = True
			else:
				mrnaAr.append( tmp )
		# end else
	# end for line
			
	# end for line
	inFile.close()
	outFile.close()
	print( 'Done' )

def handleNoncoding_celegans( inAr ):
	gene = inAr[0]
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		outStr += handleOtherCoding_celegans( mrna )
	return outStr

def handleGene_celegans( inAr, strand ):
	#print( inAr )
	
	gene = inAr[0]
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)

	for mrna in mrnas:
		# transposon coding
		if 'transposon_protein_coding' in gene['notes']:
			outStr += handleOtherCoding_celegans( mrna, replaceStr = 'transposon_transcript' )
		# other types of genes
		elif mrna[0]['feature'] != 'mRNA':
			outStr += handleOtherCoding_celegans( mrna )
		
		# regular genes
		elif strand == '+':
			outStr += handlePositive_celegans( mrna )	
		else:
			outStr += handleNegative_celegans( mrna )
	# end for
	return outStr

def handleOtherCoding_celegans( inAr, replaceStr = None ):
	mRNA = inAr[0]
	exons = inAr[1:]
	if replaceStr != None:
		mRNA['feature'] = replaceStr
	outStr = outputGene( mRNA )
	for feat in exons:
		outStr += outputGene( feat )
	return outStr
		
def handlePositive_celegans( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [feat]
	#print( mrnaID)	
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetPositive_celegans( mrnaID, UTR5, CDS, UTR3 )
	
	return outStr
	
def createGeneSetPositive_celegans( mrnaID, UTR5, CDS, UTR3 ):
	
	outStr = ''
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = i+1
		notes = 'ID=five_prime_utr:'+mrnaID+'_'+str(j) + ';'
		x = UTR5[i]
		#x['notes'] += ';' + notes
		x['notes'] = notes + x['notes']
		outStr += outputGene( x )
	
	# loop through CDS
	for i in range(len(CDS)):
		j = i+1
		x = CDS[i]
		x['notes'] = updateCDS_celegans( x['notes'], j )
		outStr += outputGene( x )
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = i+1
		notes = 'ID=three_prime_utr:'+mrnaID+'_'+str(j) + ';'
		x = UTR3[i]
		#x['notes'] += ';' + notes
		x['notes'] = notes + x['notes']
		outStr += outputGene( x )
	return outStr

def updateCDS_celegans( notes, count ):
	jndex = notes.find( 'ID=' )
	adJndex = jndex + len( 'ID=' )
	endJndex = notes[adJndex:].find( ';' )
	if endJndex == -1:
		return notes + '_{:d}'.format( count )
	else:
		newJndex = endJndex+adJndex
		return notes[:newJndex] + '_{:d}'.format( count ) + notes[newJndex:]
	
def handleNegative_celegans( inAr ):
	mRNA = inAr[0]
	# clean up notes of mrna
	outStr = outputGene( mRNA )
	mrnaID = getID( mRNA['notes'], 'ID=' )
	
	UTR5 = []
	CDS = []
	UTR3 = []
	
	for feat in inAr[1:]:
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [ feat ]
		
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetNegative_celegans( mrnaID, UTR5, CDS, UTR3 )
	return outStr

def createGeneSetNegative_celegans( mrnaID, UTR5, CDS, UTR3 ):
	outStr = ''
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = len(UTR3) - i
		notes =  'ID=three_prime_utr:'+mrnaID+'_'+str(j) + ';'
		x = UTR3[i]
		x['notes'] = notes + x['notes']
		outStr += outputGene( x )
		
	# loop through CDS
	for i in range(len(CDS)):
		j = len(CDS) - i
		x = CDS[i]
		x['notes'] = updateCDS_celegans( x['notes'], j )
		outStr += outputGene( x )
	
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = len(UTR5) - i
		notes =  'ID=five_prime_utr:'+mrnaID+'_'+str(j) + ';'
		x = UTR5[i]
		x['notes'] = notes + x['notes']
		outStr += outputGene( x )
		
	return outStr

def fixOikopleura( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	mrnaAr = []
	strand = None
	for line in inFile:
		line = line.rstrip()
		if 'gene' in line:
			#print( line  )
			# write and reset
			if len(geneAr) > 0:
				if len(mrnaAr) > 0:
					geneAr.append( mrnaAr )
				tmpStr = handleGene_Oikopleura( geneAr, strand )
				outFile.write( tmpStr )
			geneAr = []
			mrnaAr = []
			strand = None
		elif line.startswith( '#' ) :
			continue
		lineAr = line.split( '\t' )
		if len(lineAr ) < 9:
			continue
		#(0) chrm (1) source (2) feature (3) start (4) stop (5) score
		#(6) strand (7) frame (8) notes 

		tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
		if tmp['feature'] == 'gene':
			geneAr.append( tmp )
			strand = tmp['strand']
		elif tmp['feature'] == 'mRNA':
			if len(mrnaAr) > 0:
				geneAr.append( mrnaAr )
			mrnaAr = [tmp]
		else:
			mrnaAr.append( tmp )
		# end else
	# end for line
	# handle last gene
	if len(geneAr) > 0:
		if len(mrnaAr) > 0:
			geneAr.append( mrnaAr )
			tmpStr = handleGene_Oikopleura( geneAr, strand )
			outFile.write( tmpStr )
	inFile.close()
	outFile.close()
	print( 'Done' )
	
def handleGene_Oikopleura( inAr, strand ):
	gene = inAr[0]
	# fix gene notes
	geneId = searchNotes( gene['notes'], 'ID=' )
	tmp = searchNotes( gene['notes'], 'Name=' ).replace(' [Oikopleura dioica]', '')
	descAr = tmp.split(':')
	desc = (descAr[1] if len(descAr) == 2 else descAr[0])
	gene['notes'] = 'ID={:s};Name={:s};Description={:s}'.format(geneId, geneId, desc)
	mrnas = inAr[1:]
	outStr = '####\n'
	outStr += outputGene(gene)
	for mrna in mrnas:
		if strand == '+':
			outStr += handlePositive_Oikopleura( mrna, geneId )	
		else:
			outStr += handleNegative_Oikopleura( mrna, geneId )
	# end for
	return outStr

def handlePositive_Oikopleura( inAr, geneId ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mrnaID = geneId + '_T1'
	mRNA['notes'] = 'ID={:s};Parent={:s}'.format( mrnaID, geneId )
	outStr = outputGene( mRNA )
	UTR5 = []
	CDS = []
	UTR3 = []
	fivePrime = True
	for feat in inAr[1:]:
		if feat['feature'] == 'CDS':
			CDS += [ feat ]
			fivePrime = False
		if feat['feature'] == 'UTR':
			if fivePrime:
				feat['feature'] = 'five_prime_UTR'
				UTR5 += [ feat ]
			else:
				feat['feature'] = 'three_prime_UTR'
				UTR3 += [ feat ]
	#print( mrnaID)	
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetPositive_Oikopleura( mrnaID, UTR5, CDS, UTR3 )
	return outStr

def createGeneSetPositive_Oikopleura( mrnaID, UTR5, CDS, UTR3 ):
	
	outStr = ''
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = i+1
		x = UTR5[i]
		# mrnaID+'_five.prime.UTR'+str(j)
		x['notes'] = 'ID={:s}_five.prime.UTR{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	
	# loop through CDS
	for i in range(len(CDS)):
		j = i+1
		x = CDS[i]
		x['notes'] = 'ID={:s}_CDS{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = i+1
		x = UTR3[i]
		x['notes'] = 'ID={:s}_three.prime.UTR{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	return outStr

def handleNegative_Oikopleura( inAr, geneId ):
	mRNA = inAr[0]
	# clean up notes of mrna
	mrnaID = geneId + 'T1'
	mRNA['notes'] = 'ID={:s};Parent={:s}'.format( mrnaID, geneId )
	outStr = outputGene( mRNA )
	UTR5 = []
	CDS = []
	UTR3 = []
	threePrime = True
	for feat in inAr[1:]:
		if feat['feature'] == 'CDS':
			CDS += [ feat ]
			threePrime = False
		if feat['feature'] == 'UTR':
			if threePrime:
				feat['feature'] = 'three_prime_UTR'
				UTR3 += [ feat ]
			else:
				feat['feature'] = 'five_prime_UTR'
				UTR5 += [ feat ]
	#print( mrnaID)	
	# create gene/mRNA and add additional attribute information
	outStr += createGeneSetNegative_Oikopleura( mrnaID, UTR5, CDS, UTR3 )
	return outStr

def createGeneSetNegative_Oikopleura( mrnaID, UTR5, CDS, UTR3 ):
	
	outStr = ''
	# loop through UTR3
	for i in range(len(UTR3) ):
		j = len(UTR3) - i
		x = UTR3[i]
		x['notes'] = 'ID={:s}_three.prime.UTR{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	
	# loop through CDS
	for i in range(len(CDS)):
		j = len(CDS) - i
		x = CDS[i]
		x['notes'] = 'ID={:s}_CDS{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	
	# loop through UTR5
	for i in range(len(UTR5) ):
		j = len(UTR5) - i
		x = UTR5[i]
		x['notes'] = 'ID={:s}_five.prime.UTR{:d};Parent={:s}'.format( mrnaID, j, mrnaID )
		outStr += outputGene( x )
	return outStr
	
def outputGene( inDict ):
	#print( inDict )
	tmp = [ inDict['chrm'], inDict['source'], inDict['feature'], str(inDict['start']), str(inDict['end']), inDict['score'], inDict['strand'], inDict['frame'], inDict['notes'] ]
	return '{:s}\n'.format( '\t'.join( tmp ) )
