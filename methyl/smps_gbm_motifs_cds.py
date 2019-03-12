import sys, math, glob, multiprocessing, subprocess, os, functools

# Usage: python3.4 smps_gbm_motifs.py <smp_file> <fasta_file> <gff_file>

CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def readFasta( fastaFileStr ):
	''' 
		we will read in the whole fasta because it shouldn't be too big
		read into dictionary where chrm is the key and and the value is the
		sequence
		NOTE: we are making the sequence 1-based indexed and automatically capitalizing
		all of the bases
		returns the created dictionary
	'''
	fastaFile = open( fastaFileStr, 'r' )
	fastaDict = {}
	chrm = None
	seq = [' ']
	for line in fastaFile:
		line = line.rstrip()
		# chromosome headers
		if line.startswith('>'):
			# need to write old sequence
			if seq != [' '] and chrm != None:
				fastaDict[chrm] = seq
				seq = [' ']
			lineAr = line.split(' ')
			chrm = lineAr[0].replace('>', '')
		
		# sequence
		else:
			seqL = list(line.upper())
			seq += seqL
			
	# handle last chrm read
	if chrm not in fastaDict:
		fastaDict[chrm] = seq
	fastaFile.close()
	return fastaDict

def readGFF( gffFileStr ):
	
	gffFile = open( gffFileStr, 'r' )
	gffDict = {}
	
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		chrm = lineAr[0]
		start = int(lineAr[3])
		end = int(lineAr[4])
		
		if lineAr[2] == 'CDS':
			# check primary transcript
			if isPrimary( lineAr[8] ) == False:
				continue
			#check chrm
			if gffDict.get( chrm ) == None:
				gffDict[ chrm ] = []
			
			# add gene
			gffDict[chrm] += [ (start,end) ]
	# end for
	gffFile.close()
	return gffDict

def isPrimary (notesStr):
	search = "Parent="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(',')
	name = notesStr[adIndex:endIndex+adIndex]
	return name.endswith( '.1' )

def calculateExpected( gffDict, fastaDict ):
	
	nucAr = ['A','C','G','T']
	tetraDict = {}
	# get all the possible tetramers
	for n1 in nucAr:
		for n2 in nucAr:
			tetramer = n1 + 'CG' + n2
			tetraDict[tetramer] = [0,0,0]
	
	# iterate through chromosomes
	for chrm in CHRMLIST:
		fastaAr = fastaDict[chrm]
		# iterate through genes
		for gene in gffDict[chrm]:
			start = gene[0]
			end = gene[1]
			# iterate through sequence
			for i in range( start, end-2 ):
				try:

					if fastaAr[i+1] == 'C' and fastaAr[i+2] == 'G':
						# positive strand
						tetra = ''.join( fastaAr[i:i+4] )
						tetraDict[ tetra ][0] += 1
						# negative strand
						tmp = reverseComplement(fastaAr[i:i+4])
						tetra = ''.join( tmp )
						tetraDict[ tetra ][0] += 1
				except KeyError:
					pass
	return tetraDict

def reverseComplement( seq ):
	nucAr = ['A','C','G','T']
	revAr = ['T','G','C','A']
	newSeq = []
	
	# iterate over bases
	for i in range( len(seq) ):
		base = seq[ len(seq)-1-i ]
		ind = nucAr.index( base )
		newSeq += [ revAr[ind] ]
	return newSeq

def readSMP( smpFileStr ):
	'''
		reads the file containing CG-SMPs
		tab-delimited file with 2 columns: chromosome and position
	'''

	smpFile = open( smpFileStr, 'r' )
	smpDict = {}
	
	for line in smpFile:
		if line.startswith( '#' ) or line == '\n':
			continue
		lineAr = line.rstrip().split( '\t' )
		
		chrm = lineAr[0]
		pos = int( lineAr[1] )
		# chrm not in dict
		if smpDict.get( chrm ) == None:
			smpDict[ chrm ] = []
		smpDict[ chrm ] += [pos]
	smpFile.close()
	return smpDict

def calculateObserved( smpDict, fastaDict, tetraDict ):
	
	# iterate through chromosomes
	for chrm in smpDict.keys():
		# iterate through smp positions
		seq = fastaDict[chrm]
		for pos in smpDict[chrm]:
			b = seq[pos]
			if b == 'C':
				# n1 C G n2
				n1 = seq[pos-1]
				n2 = seq[pos+2]
			elif b == 'G':
				# rc(n2) G C rc(n1)
				n1 = baseComplement(seq[pos+1])
				n2 = baseComplement(seq[pos-2])
			else:
				print( 'ERROR: {:s} {:d} CG-SMP does not align with proper cytosine...ignoring'.format( chrm, pos ) )
				continue
			# determine tetramer
			tetra = n1 + 'C' + 'G' + n2
			# add to tetraDict
			tetraDict[tetra][1] += 1
		# end for ps
	# end for chrm
	return tetraDict

def baseComplement( base ):
	nucAr = ['A','C','G','T']
	revAr = ['T','G','C','A']
	ind = nucAr.index( base )
	return revAr[ind]

def calculateCorrected( tetraDict ):
	'''
		fills in the last array position for tetraDict by calculating the
		corrected frequency
		corrected frequency = observed / expected
	'''
	
	for tetra in tetraDict.keys():
		tetraDict[tetra][2] = float( tetraDict[tetra][1] ) / tetraDict[tetra][0]
	return tetraDict

def writeOutput( outFileStr, tetraDict ):
	'''
		writes output to outFileStr
		tab-delimited file with 4 columns: tetramer, number expected,
		number observed, adjusted
	'''
	
	outFile = open( outFileStr, 'w' )
	# header
	outFile.write( '#NCGN\texpected\tobserved\tadjusted\n' )
	# iterate through tetramers
	for tetramer in sorted( tetraDict.keys() ):
		vals = tetraDict[tetramer]
		outStr = '{:s}\t{:d}\t{:d}\t{:.6f}\n'.format( tetramer, vals[0], vals[1], vals[2] )
		outFile.write( outStr )
	outFile.close()

def processInputs( smpFileStr, fastaFileStr, gffFileStr ):
	'''
		main function
		read FASTA, calculate expected tetramer content, find observed tetramer
		number, adjust, then write to outfile
	'''
	
	outFileStr = smpFileStr[:(smpFileStr.rfind('.'))] + '_ncgn_cds.tsv'
	print( 'Reading FASTA...' )
	fastaDict = readFasta( fastaFileStr )
	print( 'Reading GFF and counting number of tetramers in genes...' )
	gffDict = readGFF( gffFileStr )
	print( 'Calculating expected number of counts for tetramers...' )
	tetraDict = calculateExpected( gffDict, fastaDict )
	print( 'Reading CG-SMPs file...' )
	smpDict = readSMP( smpFileStr )
	print( 'Calculating observed counts for tetramers...' )
	tetraDict = calculateObserved( smpDict, fastaDict, tetraDict )
	tetraDict = calculateCorrected( tetraDict )
	print( 'Writing output to {:s}...'.format( outFileStr ) )
	writeOutput( outFileStr, tetraDict )
	
if __name__ == "__main__":
	if len(sys.argv) != 4 :
		print ("Usage: python3.4 smps_gbm_motifs.py <smp_file> <fasta_file> <gff_file>")
	else:
		smpFileStr = sys.argv[1]
		fastaFileStr = sys.argv[2]
		gffFileStr = sys.argv[3]
		processInputs( smpFileStr, fastaFileStr, gffFileStr )
