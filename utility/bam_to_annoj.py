import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 bam_to_annoj.py <bam/sam_file> [bam/sam_file]*

'''
+----------+--------+---------+---------+--------------------------------------------+-----------+----+
| assembly | strand | start   | end     | sequenceA                                  | sequenceB | id |
+----------+--------+---------+---------+--------------------------------------------+-----------+----+
| 6        | -      | 3000096 | 3000137 | AAAAGGGCTGAGTTTCCTGTATTCGAGGTCTTCCTGATAAGT | NULL      |  1 | 
| 6        | +      | 3000130 | 3000171 | CCCTTTTATAACCCACAAAGAGTGGTTGGAGACGTGTCAGCT | NULL      |  2 | 
| 6        | +      | 3000130 | 3000171 | CCCTTTTATAACCCACAAAGAGTGGTTGGAGACGTGTCAGCT | NULL      |  3 | 
| 6        | +      | 3000204 | 3000245 | CTGCCACGGGACAGGCAGGGCACAAGGCATGGAAAAATACCC | NULL      |  4 | 
| 6        | -      | 3000308 | 3000349 | CAGTCAACCCTCTGTGTGAGGAGCCCCCCTTGCAATCGCCAT | NULL      |  5 | 
+----------+--------+---------+---------+--------------------------------------------+-----------+——+

SEPARATE FILE FOR EACH CHROMOSOME
'''

def checkBamFile( bamFileStr ):
	
	# if bam, convert to sam
	# samtools must be part of path
	if bamFileStr.endswith( '.bam' ):
		outSamFileStr = bamFileStr.replace('.bam','.sam' )
		bamCommand = "samtools view -h {:s} > {:s}".format( bamFileStr, outSamFileStr )
		subprocess.call( bamCommand, shell=True )
		return outSamFileStr
	else:
		return False

def readSam( samFileStr ):

	samFile = open( samFileStr, 'r' )
	outPre = determineOutFilePrefix( samFileStr )
	currChrm = None
	currFile = None
	
	# look through sam file
	for line in samFile:
		# headers
		line=line.rstrip()
		if line.startswith('@'):
			continue
		lineAr = line.split('\t')
		# (0) name (1) flag (2) chromosome (3) 1-based position (4) mapping 
		# quality (5) CIGAR string (6) chromosome of mate (7) position of mate 
		# (8) sign-observed template length (9) sequence (10) base quality +33 
		# (11-) optional fields
		
		readStart = int( lineAr[3] )
		seq = lineAr[9]
		readEnd = readStart + len( seq )
		chrm = lineAr[2]
		if chrm.startswith('scaffold'):
			continue
		
		# decode flag
		# [0] read is 2nd mate [1] read is 1st mate [2] mate on reverse strand
		# [3] read on reverse strand [4] mate unmapped [5] read unmapped [6] read mapped
		# in proper pair [7] read is paired
		flag = '{:08b}'.format( int( lineAr[1] ) )
		strand = None
		# mapped and mapped to reverse strand
		if flag[5] == '0' and flag[3] == '1':
			strand = '-'
		# mapped and not mapped to reverse strand, so mapped to forward strand
		elif flag[5] == '0':
			strand = '+'
		
		# check curr chrm
		if chrm == '*':
			continue
		elif currChrm != chrm:
			# close previous
			if currFile != None:
				currFile.close()
			# change currChrm
			currChrm = chrm
			# create curFile
			currFileStr = '{:s}_annoj_{:s}.txt'.format( outPre, chrm )
			currFile = open( currFileStr, 'w' )
		# write line
		# assembly strand start end sequence NULL
		#print( currChrm, strand, readStart, readEnd, seq )
		try:
			outStr = '{:s}\t{:s}\t{:d}\t{:d}\t{:s}\tNULL\n'.format( currChrm, strand, readStart, readEnd, seq )
		except TypeError:
			#print( currChrm, strand, readStart, readEnd, seq )
			pass
		currFile.write( outStr )
	
	# close files
	currFile.close()
	samFile.close()
		
def determineOutFilePrefix( samFileStr ):
	rind = samFileStr.rfind( '.' )
	return samFileStr[:rind]

def mainFunction( bamFileStrAr ):
	
	for bamFileStr in bamFileStrAr:
		print( 'Working on {:s}...'.format( bamFileStr ) )
		
		samFileStr = checkBamFile( bamFileStr )
		if samFileStr == False:
			readSam( bamFileStr )
		else:
			readSam( samFileStr )
			os.remove( samFileStr )
	print( 'Done' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3.4 bam_to_annoj.py <bam/sam_file> [bam/sam_file]")
	else:
		bamFileStrAr = []
		for i in range( 1, len(sys.argv) ):
			bamFileStrAr += [sys.argv[i]]
		mainFunction( bamFileStrAr )
