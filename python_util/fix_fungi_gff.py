import math, copy, functools, re

def reformatInitial( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	currGene = None
	for line in inFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) source (2) feature (3) start (4) end (5) score
		# (6) strand (7) phase (8) attributes
		# get transcript name
		gName = gtfSearch( lineAr[8], 'name' )
		# same gene
		if gName == currGene:
			outFile.write( line )
		# different gene
		else:
			outFile.write( '###\n' )
			outFile.write( line )
			currGene = gName
	# end for line
	outFile.write( '###\n' )
	inFile.close()
	outFile.close()

def quickReformat( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	currGene = None
	for line in inFile:
		if line.startswith('#'):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chrm (1) source (2) feature (3) start (4) end (5) score
		# (6) strand (7) phase (8) attributes
		if lineAr[2] == 'gene':
			outFile.write( '###\n' )
			outFile.write( line )
		else:
			outFile.write( line )
	# end for line
	outFile.write( '###\n' )
	inFile.close()
	outFile.close()

def gtfSearch( inStr, searchStr ):
	spAr = inStr.split(';')
	for sp in spAr:
		# name "gene name"
		# transcriptId idNum
		#print(sp)
		spAr2 = sp.lstrip().rstrip().split(' ')
		if spAr2[0] == searchStr:
			return spAr2[1].replace('"', '')
	# end for sp
	return False

def reformatJGI( inFileStr, outFileStr, genePre, speciesPortal ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	for line in inFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		# write and reset
		if line == '###':
			if len(geneAr) > 0:
				tmpStr = handleGeneJGI( geneAr, genePre, speciesPortal )
				outFile.write( tmpStr )
			geneAr = []
		elif line.startswith('#') or len(lineAr) < 9:
			continue
		else:
			lineAr = line.rstrip().split('\t')
			# (0) chrm (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) phase (8) attributes
			# get transcript name
			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			geneAr += [tmp]
	# end for line
	inFile.close()
	outFile.close()

def handleGeneJGI( geneAr, genePre, speciesPortal ):
	# basic info we need
	strand = None
	geneName = None
	geneId = None
	proteinId = None
	for feat in geneAr:
		if strand == None:
			strand = feat['strand']
		gName = gtfSearch(feat['notes'], 'name' )
		if geneName == None and gName:
			geneName = gName
		gId = gtfSearch(feat['notes'], 'transcriptId' )
		if geneId == None and gId:
			geneId = gId
		pId = gtfSearch(feat['notes'], 'proteinId' )
		if proteinId == None and pId:
			proteinId = pId
	# process the info we got based on strand
	#print( geneAr, geneName, geneId)
	if strand == '+':
		return handlePlusGene( geneAr, speciesPortal, genePre, geneName, geneId, proteinId )
	else:
		return handleMinusGene( geneAr, speciesPortal, genePre, geneName, geneId, proteinId )

def handlePlusGene( inGeneAr, speciesPortal, genePre, geneName, geneId, proteinId ):
	
	UTR5 = []
	CDS = []
	UTR3 = []
	cnt = 0
	fStart = -1
	fEnd = -1
	start = -1
	end = -1
	# sort geneAr
	geneAr = sorted( inGeneAr, key=functools.cmp_to_key(customSortPlus) )
	
	# get gene info
	geneStart = math.inf
	geneEnd = -math.inf
	geneChrm = geneAr[0]['chrm']
	geneSource = geneAr[0]['source']
	
	for feat in geneAr:
		# update total gene info
		if feat['start'] < geneStart:
			geneStart = feat['start']
		if feat['end'] > geneEnd:
			geneEnd = feat['end']
		# start codon -> if UTR/exon, adjust
		if feat['feature'] == 'start_codon':
			start = feat['start'] - 1
			cnt += 1
			if len(UTR5) > 0:
				tmp = UTR5[-1]
				fEnd = tmp['end']
				tmp['end'] = start
		# stop codon -> if we have CDS
		elif feat['feature'] == 'stop_codon':
			end = feat['end'] + 1
			cnt += 1
			tmp = CDS[-1]
			# single exon
			if tmp['feature'] != 'exon':
				tmp = copy.copy(UTR5[-1])
				tmp['end'] = fEnd
			# multiple exons
			else:
				CDS = CDS[:-1]
			tmp['start'] = end
			if tmp['start'] < tmp['end']:
				UTR3 = [tmp]
		# UTR5 exons
		elif feat['feature'] == 'exon' and cnt == 0:
			UTR5 += [feat]
		# CDS exons
		elif feat['feature'] == 'exon' and cnt == 1:
			CDS += [feat]
		# CDS features
		elif feat['feature'] == 'CDS':
			CDS += [feat]
		# UTR3 exons
		elif feat['feature'] == 'exon' and cnt == 2:
			UTR3 += [feat]
	# end for feat
	
	# clean up CDS feature list
	newCDS = []
	for x in CDS:
		if x['feature'] == 'CDS':
			newCDS += [x]
	
	geneBasics = {'start': geneStart, 'end': geneEnd, 'chrm': geneChrm, 'source': geneSource, 'score':'.', 'strand': '+', 'frame': '.' }
	
	return createGeneSetPlus( speciesPortal, genePre, geneName, geneId, proteinId, geneBasics, UTR5, newCDS, UTR3 )
					
def handleMinusGene( inGeneAr, speciesPortal, genePre, geneName, geneId, proteinId ):
	UTR5 = []
	CDS = []
	UTR3 = []
	cnt = 0
	fStart = -1
	fEnd = -1
	start = -1
	end = -1
	# sort geneAr
	geneAr = sorted( inGeneAr, key=functools.cmp_to_key(customSortMinus) )
	
	# get gene info
	geneStart = math.inf
	geneEnd = -math.inf
	geneChrm = geneAr[0]['chrm']
	geneSource = geneAr[0]['source']
	#print( '####', genePre, geneName, geneId )
	#print( geneAr )
	for feat in geneAr:
		# update total gene info
		if feat['start'] < geneStart:
			geneStart = feat['start']
		if feat['end'] > geneEnd:
			geneEnd = feat['end']
		# stop codon -> if UTR/exon, adjust
		if feat['feature'] == 'stop_codon':
			start = feat['start'] - 1
			cnt += 1
			if len(UTR3) > 0:
				tmp = UTR3[-1]
				fEnd = tmp['end']
				tmp['end'] = start
		# start codon
		elif feat['feature'] == 'start_codon':
			end = feat['end'] + 1
			cnt += 1
			tmp = CDS[-1]
			# single exon for gene
			if tmp['feature'] != 'exon':
				tmp = copy.copy(UTR3[-1])
				tmp['end'] = fEnd
			# start in middle/end of exon
			else:
				CDS = CDS[:-1]
			tmp['start'] = end
			if tmp['start'] < tmp['end']:
				UTR5 = [tmp]
		# UTR3 exons
		elif feat['feature'] == 'exon' and cnt == 0:
			UTR3 += [feat]
		# CDS exons
		elif feat['feature'] == 'exon' and cnt == 1:
			CDS += [feat]
		# CDS features
		elif feat['feature'] == 'CDS':
			CDS += [feat]
		# UTR5 exons
		elif feat['feature'] == 'exon' and cnt == 2:
			UTR5 += [feat]
	# end for feat
	
	# clean up CDS feature list
	newCDS = []
	for x in CDS:
		if x['feature'] == 'CDS':
			newCDS += [x]
	
	geneBasics = {'start': geneStart, 'end': geneEnd, 'chrm': geneChrm, 'source': geneSource, 'score':'.', 'strand': '-', 'frame': '.' }
	#print( UTR5, newCDS, UTR3)
	out = createGeneSetMinus( speciesPortal, genePre, geneName, geneId, proteinId, geneBasics, UTR5, newCDS, UTR3 )
	#print( out )
	return out

def customSortPlus( featA, featB ):
	featList = ['start_codon', 'CDS', 'exon', 'stop_codon']
	if featA['start'] > featB['start']:
		return 1
	elif featA['start'] == featB['start']:
		try:
			aIndex = featList.index(featA['feature'])
			bIndex = featList.index(featB['feature'])
			if aIndex == bIndex:
				return 0
			elif aIndex > bIndex:
				return 1
			else:
				return -1
		except ValueError:
			print('WARNING: feature {:s} or {:s} not found'.format( featA['feature'], featB['feature']) )
	else:
		return -1

def createGeneSetPlus( speciesPortal, genePre, geneName, geneId, proteinId, geneBasics, UTR5, CDS, UTR3 ):
	outStr = '###\n'
	# ID=gene_286;Name=jgi.p|Aspfl1|25780;portal_id=Aspfl1
	formatGeneId = 'gene_{:s}'.format( geneId )
	formatGeneName = '{:s}.{:s}'.format( genePre, geneId )
	formatmRNAId = 'mRNA_{:s}'.format( geneId )
	
	# create gene
	geneFeat = copy.copy( geneBasics )
	geneFeat['feature'] = 'gene'
	geneFeat['notes'] = 'ID={:s};Name={:s};Portal_id={:s};Transcript_id={:s};Protein_id={:s}'.format( formatGeneId, formatGeneName, speciesPortal, geneId, proteinId )
	if geneName != '' and geneName != None:
		geneFeat['notes'] += ';Description={:s}'.format( geneName )
	
	outStr += outputFeature( geneFeat )
	
	# create mRNA
	mrnaFeat = copy.copy( geneBasics )
	mrnaFeat['feature'] = 'mRNA'
	mrnaFeat['notes'] = 'ID={:s};Name={:s}.t;Parent={:s};Transcript_id={:s};Protein_id={:s}'.format( formatmRNAId, formatGeneName, formatGeneId, geneId, proteinId )
	outStr += outputFeature( mrnaFeat )
	
	# loop through UTR5
	for i in range(len(UTR5)):
		outStr += outputCountedFeature( UTR5[i], formatmRNAId, i+1, 'five_prime_UTR', 'UTR5')
	
	# loop through CDS
	for i in range(len(CDS)):
		outStr += outputCountedFeature( CDS[i], formatmRNAId, i+1, 'CDS', 'CDS')
	
	# loop through UTR3
	for i in range(len(UTR3)):
		outStr += outputCountedFeature( UTR3[i], formatmRNAId, i+1, 'three_prime_UTR', 'UTR3')
	
	return outStr

def customSortMinus( featA, featB ):
	featList = ['stop_codon', 'CDS', 'exon', 'start_codon']
	if featA['start'] > featB['start']:
		return 1
	elif featA['start'] == featB['start']:
		try:
			aIndex = featList.index(featA['feature'])
			bIndex = featList.index(featB['feature'])
			if aIndex == bIndex:
				return 0
			elif aIndex > bIndex:
				return 1
			else:
				return -1
		except ValueError:
			print('WARNING: feature {:s} or {:s} not found'.format( featA['feature'], featB['feature']) )
	else:
		return -1

def createGeneSetMinus( speciesPortal, genePre, geneName, geneId, proteinId, geneBasics, UTR5, CDS, UTR3 ):
	# ID=gene_286;Name=jgi.p|Aspfl1|25780;portal_id=Aspfl1
	outStr = '###\n'
	formatGeneId = 'gene_{:s}'.format( geneId )
	formatGeneName = '{:s}.{:s}'.format( genePre, geneId )
	formatmRNAId = 'mRNA_{:s}'.format( geneId )
	
	# create gene
	geneFeat = copy.copy( geneBasics )
	geneFeat['feature'] = 'gene'
	geneFeat['notes'] = 'ID={:s};Name={:s};Portal_id={:s};Transcript_id={:s};Protein_id={:s}'.format( formatGeneId, formatGeneName, speciesPortal, geneId, proteinId )
	if geneName != '' and geneName != None:
		geneFeat['notes'] += ';Description={:s}'.format( geneName )
	outStr += outputFeature( geneFeat )
	
	# create mRNA
	mrnaFeat = copy.copy( geneBasics )
	mrnaFeat['feature'] = 'mRNA'
	mrnaFeat['notes'] = 'ID={:s};Name={:s}.t;Parent={:s};Transcript_id={:s};Protein_id={:s}'.format( formatmRNAId, formatGeneName, formatGeneId, geneId, proteinId  )
	outStr += outputFeature( mrnaFeat )
	
	# loop through UTR3
	for i in range(len(UTR3)):
		outStr += outputCountedFeature( UTR3[i], formatmRNAId, len(UTR3) - i, 'three_prime_UTR', 'UTR3')
	
	# loop through CDS
	for i in range(len(CDS)):
		outStr += outputCountedFeature( CDS[i], formatmRNAId, len(CDS) - i, 'CDS', 'CDS')
	
	# loop through UTR5
	for i in range(len(UTR5)):
		outStr += outputCountedFeature( UTR5[i], formatmRNAId, len(UTR5) - i, 'five_prime_UTR', 'UTR5')
	
	return outStr


def quickFormatJGI( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	for line in inFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		if line == '###':
			if len(geneAr) > 0:
				tmpStr = handleQuickGeneJGI( geneAr )
				outFile.write( tmpStr )
			geneAr = []
		elif line.startswith('#') or len(lineAr) < 9:
			continue
		else:
			lineAr = line.rstrip().split('\t')
			# (0) chrm (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) phase (8) attributes
			# get transcript name
			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			geneAr += [tmp]
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

def handleQuickGeneJGI( geneAr ):
	outStr = '###\n'
	gene = geneAr[0]
	strand = gene['strand']
	# update gene notes
	geneId = searchNotes( gene['notes'], 'ID=' )
	geneName = searchNotes( gene['notes'], 'Name=' )
	genePortal  = searchNotes( gene['notes'], 'portal_id=' )
	transcriptId = searchNotes( gene['notes'], 'transcriptId=' )
	proteinId = searchNotes( gene['notes'], 'proteinId=' )
	formatGeneName = geneName.replace('jgi.p|', '').replace('|', '.')
	
	gene['notes'] = 'ID={:s};Name={:s};Portal_id={:s};Transcript_id={:s};Protein_id={:s}'.format( geneId, formatGeneName, genePortal, transcriptId, proteinId )
	outStr += outputFeature( gene )
	
	mRNA = geneAr[1]
	mRNAId = searchNotes( mRNA['notes'], 'ID=' )
	mRNA['notes'] = 'ID={:s};Name={:s}.t;Parent={:s};Transcript_id={:s};Protein_id={:s}'.format( mRNAId, formatGeneName, geneId, transcriptId, proteinId )
	outStr += outputFeature( mRNA )
	
	# handle other features by strand
	if strand == '+':
		outStr += handleQuickGenePlus( geneAr[2:], mRNAId )
	else:
		outStr += handleQuickGeneMinus( geneAr[2:], mRNAId )
	return outStr
	
def handleQuickGenePlus( geneAr, mRNAId ):
	UTR5 = []
	CDS = []
	UTR3 = []
	for feat in geneAr:
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [feat]
	
	# output
	outStr = ''
	for i in range(len(UTR5)):
		outStr += outputCountedFeature( UTR5[i], mRNAId, i+1, 'five_prime_UTR', 'UTR5' )
	for i in range(len(CDS)):
		outStr += outputCountedFeature( CDS[i], mRNAId, i+1, 'CDS', 'CDS' )
	for i in range(len(UTR3)):
		outStr += outputCountedFeature( UTR3[i], mRNAId, i+1, 'three_prime_UTR', 'UTR3' )
	
	return outStr

def handleQuickGeneMinus( geneAr, mRNAId ):
	UTR5 = []
	CDS = []
	UTR3 = []
	for feat in geneAr:
		if feat['feature'] == 'five_prime_UTR':
			UTR5 += [ feat ]
		elif feat['feature'] == 'CDS':
			CDS += [ feat ]
		elif feat['feature'] == 'three_prime_UTR':
			UTR3 += [feat]
	
	# output
	outStr = ''
	for i in range(len(UTR3)):
		outStr += outputCountedFeature( UTR3[i], mRNAId, len(UTR3) - i, 'three_prime_UTR', 'UTR3' )
	for i in range(len(CDS)):
		outStr += outputCountedFeature( CDS[i], mRNAId, len(CDS) - i, 'CDS', 'CDS' )
	for i in range(len(UTR5)):
		outStr += outputCountedFeature( UTR5[i], mRNAId, len(UTR5) - i, 'five_prime_UTR', 'UTR5' )
	
	return outStr

def addCodons( inFileStr, outFileStr ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	geneAr = []
	for line in inFile:
		line = line.rstrip()
		lineAr = line.split('\t')
		# write and reset
		if line == '###':
			if len(geneAr) > 0:
				tmpStr = handleGeneCodon( geneAr )
				outFile.write( tmpStr )
			geneAr = []
		elif line.startswith('#') or len(lineAr) < 9:
			continue
		else:
			lineAr = line.rstrip().split('\t')
			# (0) chrm (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) phase (8) attributes
			# get transcript name
			tmp = {'chrm':lineAr[0], 'source':lineAr[1], 'feature':lineAr[2], 'start': int(lineAr[3]), 'end': int(lineAr[4]), 'score':lineAr[5], 'strand':lineAr[6], 'frame':lineAr[7], 'notes': lineAr[8] }
			geneAr += [tmp]
	# end for line
	inFile.close()
	outFile.close()

def handleGeneCodon( geneAr ):
	strand = geneAr[0]['strand']
	gName = gtfSearch(geneAr[0]['notes'], 'name' )
	
	# need to add stop codon
	if strand == '+':
		return handlePlusGeneCodon( geneAr, gName )
	else:
		return handleMinusGeneCodon( geneAr, gName )
	return outStr

def handlePlusGeneCodon( inGeneAr, geneId ):
	# sort
	geneAr = sorted( inGeneAr, key=functools.cmp_to_key(customSortPlus) )
	
	cdsAr = []
	exonAr = []
	offset = 0
	
	# split by CDS and exon
	for feat in geneAr:
		if feat['feature'] == 'CDS':
			cdsAr += [feat]
		elif feat['feature'] == 'exon':
			exonAr += [feat]
	
	# add start codon
	tmp = copy.copy(cdsAr[0])
	tmp['end'] = tmp['start'] + 2
	tmp['feature'] = 'start_codon'
	tmp['notes'] = 'name "{:s}"'.format( geneId )
	geneAr += [tmp]
	
	# if exon CDS lengths are different, there's utrs to account for -> account for 5'
	if len(cdsAr) < len(exonAr):
		for i in range(len(exonAr)):
			for j in range(len(cdsAr)):
				if exonAr[i]['start'] == cdsAr[j]['start'] or exonAr[i]['end'] == cdsAr[j]['end']:
					offset = i-j
					break
			# end for j
		# end for i
		#print( 'offset:', offset ,'CDS:', len(cdsAr), 'Exon:', len(exonAr), 'Gene:', cdsAr[0]['notes'] )
	elif len(cdsAr) > len(exonAr):
		print( 'ERROR: CDS: ', len(cdsAr), 'Exon: ', len(exonAr), 'Gene:', cdsAr[0]['notes'] )
		exit()
	
	# otherwise, loop through pairs backwards
	endPos = None
	maxCoord = -1
	for i in range(len(cdsAr)-1, -1, -1):
		tmpCDS = cdsAr[i]
		tmpExon = exonAr[i+offset]
		maxCoord = max(maxCoord, tmpCDS['end'], tmpExon['end'] )
		# found end
		if tmpExon['end'] > tmpCDS['end']:
			endPos = tmpCDS['end']
	
	# we found end with UTR
	tmp = copy.copy(exonAr[0])
	if endPos == None:
		endPos = maxCoord
	# update feature
	tmp['end'] = endPos
	tmp['start'] = endPos - 2
	tmp['feature'] = 'stop_codon'
	tmp['notes'] = 'name "{:s}"'.format( geneId )
	geneAr += [tmp]
	
	# reorder again
	orderedGeneAr = sorted( geneAr, key=functools.cmp_to_key(customSortPlus) )
	outStr = '###\n'
	for feat in orderedGeneAr:
		outStr += outputFeature( feat )
	return outStr

def handleMinusGeneCodon( inGeneAr, geneId ):
	# sort
	geneAr = sorted( inGeneAr, key=functools.cmp_to_key(customSortMinus) )
	
	cdsAr = []
	exonAr = []
	offset = 0
	
	# split by CDS and exon
	for feat in geneAr:
		if feat['feature'] == 'CDS':
			cdsAr += [feat]
		elif feat['feature'] == 'exon':
			exonAr += [feat]
	
	# add start codon
	tmp = copy.copy(cdsAr[-1])
	tmp['start'] = tmp['end'] - 2
	tmp['feature'] = 'start_codon'
	tmp['notes'] = 'name "{:s}"'.format( geneId )
	geneAr += [tmp]
	
	# if exon CDS lengths are different, there's a problem
	if len(cdsAr) < len(exonAr):
		for i in range(len(exonAr)):
			for j in range(len(cdsAr)):
				if exonAr[i]['start'] == cdsAr[j]['start'] or exonAr[i]['end'] == cdsAr[j]['end']:
					offset = i-j
					break
			# end for j
		# end for i
		#print( 'offset:', offset ,'CDS:', len(cdsAr), 'Exon:', len(exonAr), 'Gene:', cdsAr[0]['notes'] )
	elif len(cdsAr) > len(exonAr):
		print( 'ERROR: CDS: ', len(cdsAr), 'Exon: ', len(exonAr), 'Gene:', cdsAr[0]['notes'] )
		exit()
	
	# otherwise, loop through pairs backwards
	startPos = None
	minCoord = float('inf')
	for i in range(len(cdsAr)):
		tmpCDS = cdsAr[i]
		tmpExon = exonAr[i+offset]
		maxCoord = min(minCoord, tmpCDS['start'], tmpExon['start'] )
		# found start
		if tmpExon['start'] < tmpCDS['start']:
			startPos = tmpCDS['start']
	
	# we found end with UTR
	tmp = copy.copy(exonAr[0])
	if startPos == None:
		startPos = maxCoord
	# update feature
	tmp['start'] = startPos
	tmp['end'] = startPos + 2
	tmp['feature'] = 'stop_codon'
	tmp['notes'] = 'name "{:s}"'.format( geneId )
	geneAr += [tmp]
	
	orderedGeneAr = sorted( geneAr, key=functools.cmp_to_key(customSortMinus) )
	outStr = '###\n'
	for feat in orderedGeneAr:
		outStr += outputFeature( feat )
	return outStr
		
def outputCountedFeature( feat, mRNA, j, feature, fType ):
	feat['notes'] = 'ID={:s}.{:s}_{:d};Parent={:s}'.format( mRNA, fType, j, mRNA )
	feat['feature'] = feature
	return outputFeature(feat)
	
def outputFeature( inDict ):
	#print( inDict )
	tmp = [ inDict['chrm'], inDict['source'], inDict['feature'], str(inDict['start']), str(inDict['end']), inDict['score'], inDict['strand'], inDict['frame'], inDict['notes'] ]
	return '{:s}\n'.format( '\t'.join( tmp ) )

