import sys, math, os
from bioFiles import *

# Usage: python sdi_to_vcf.py <sampleName> <fasta_file> <in_file> <out_file>

def processInputs( inFileStr, outFileStr, fastaFileStr, sampleName ):
	print( 'Converting', os.path.basename(inFileStr), 'to', os.path.basename(outFileStr) )
	
	# read fasta
	fasta = FileFASTA( fastaFileStr )
	fastaDict = fasta.getFastaDict()
	
	# SDI header
	header = outHeader( sampleName )
	# read SDI
	i, s, u = readSDI( inFileStr, outFileStr, fastaDict, header )
	print( 'Ambiguous insertions:', i )
	print( 'Ambiguous SNPs:', s )
	print( 'Successful SN:', u )

def readSDI( inFileStr, outFileStr, fastaDict, header ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	outFile.write( header )
	ambigIns = 0
	ambigSnp = 0
	succ = 0
	iupac = getIUPAC()
	
	for line in inFile:
		lineAr = line.rstrip().split('\t')
		# (0) chromosome (1) position (2) length (3) reference
		# (4) consensus (5) detection score (6) HMQ coverage (depth)
		# (7) Phred score (8) HMQ consensus
		chrm, pos, length, ref, alt = lineAr[0:5]
		pos = int(pos)
		length = int( length )
		
		depth = ('0' if len(lineAr) <= 6 else lineAr[6] )
		
		# VCF format
		# (0) chrm (1) pos (2) ID (3) ref (4) alt (5) qual
		# (6) filter (7) info (8) format
		# insertion
		if length > 0:
			# ambiguous
			if alt not in ['A', 'C', 'G', 'T', 'N' ]:
				ambigIns += 1
				continue
			refChar = fastaDict[chrm][pos-1]
			outStr = '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\tPASS\tSV=INS\tGT\t1/1\n'.format( chrm, pos, refChar, refChar+alt )
			
		# deletion
		elif length < 0:
			refChar = fastaDict[chrm][pos-2]
			outStr = '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\tPASS\tSV=DEL\tGT\t1/1\n'.format( chrm, pos-1, refChar+ref, refChar )
		else: #SNP
			if ref not in ['A', 'C', 'G', 'T', 'N' ]:
				ambigSnp += 1
				continue
			tmpIupac = iupac[alt]
			# homozygous
			if len( tmpIupac ) == 1:
				outStr = '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\tPASS\tSV=SNP;DP={:s}\tGT\t1/1\n'.format( chrm, pos, ref, alt, depth )
			elif len( tmpIupac ) == 2:
				if ref not in tmpIupac:
					continue
				s = set(tmpIupac)
				s.remove(ref)
				outStr = '{:s}\t{:d}\t.\t{:s}\t{:s}\t.\tPASS\tSV=SNP;DP={:s}\tGT\t1/0\n'.format( chrm, pos, ref, list(s)[0], depth )
			else:
				ambigSnp += 1
				continue
		succ += 1
		outFile.write( outStr )
	# end for line
	inFile.close()
	outFile.close()
	return ambigIns, ambigSnp, succ	
		

def outHeader( sampleName ):

	outStr = '##fileformat=VCFv4.1\n'
	outStr += '##FORMAT=<ID=SV,Number=1,Type=String,Description="SNV Type">\n'
	outStr += '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
	outStr += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
	outStr += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{:s}\n'.format( sampleName )
	
	return outStr

def getIUPAC():
	iupac = {
	'A':set(['A']),
	'C':set(['C']),
	'G':set(['G']),
	'T':set(['T']),
	'R':set(['A','G']),
	'Y':set(['C','T']),
	'M':set(['A','C']),
	'K':set(['G','T']),
	'W':set(['A','T']),
	'S':set(['C','G']),
	'B':set(['C','G','T']),
	'D':set(['A','G','T']),
	'H':set(['A','C','T']),
	'V':set(['A','C','G']),
	'N':set(['A','C','G','T'])}
	return iupac
	
def parseInputs( argv ):
	sampleName = argv[0]
	fastaFileStr = argv[1]
	inFileStr = argv[2]
	outFileStr = argv[3]
	processInputs( inFileStr, outFileStr, fastaFileStr, sampleName )


if __name__ == "__main__":
	if len(sys.argv) < 5 :
		print ("Usage: python sdi_to_vcf.py <sampleName> <fasta_file> <in_file> <out_file>")
	else:
		parseInputs( sys.argv[1:] )
