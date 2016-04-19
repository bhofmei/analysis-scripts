import math, sys

# Usage: count_fasta.py <fast_file>

AR = ["CG", "CHG", "CHH", "CH", "TG", "THG", "THH", "TH", "C/G", "N"]
CH = ["CA", "CC", "CT"]
CHG = [ "CAG", "CCG", "CTG" ]
CHH = [ "CAA", "CAC", "CAT", "CCA", "CCC", "CCT","CTA", "CTC", "CTT"]
TH = ["TA", "TC", "TT"]
THG = [ "TAG", "TCG", "TTG" ]
THH = [ "TAA", "TAC", "TAT", "TCA", "TCC", "TCT","TTA", "TTC", "TTT"]
CG = [ "C", "G" ]

def readFasta( fastaFileStr ):
	fastaFile = open( fastaFileStr, 'r' )
	inAr = [0] * len(AR)
	seq = ""
	for line in fastaFile:
		if line.startswith( '>' ):
			# get numbers
			inAr = analyzeSequence(seq, inAr)
			seq = "";
		else:
			seq += line.rstrip();
	# get numbers again
	inAr = analyzeSequence(seq, inAr)
	fastaFile.close()
	return inAr
	
def analyzeSequence( seq, inAr ):
	# inAr = (0) cg (1) chg (2)chh (3) ch (4) tg (5) thg (6) thh
	# (7) th (8) c/g (9) n
	print( "start counting" )
	inAr[9] += len(seq)
	# loop through sequence
	for i in range(len(seq)):
		try:
			#1 character
			if seq[i] in CG:
				inAr[8] += 1
			# 2 characters
			c2 = seq[i:i+2]
			if c2 =="CG":
				inAr[0] += 1
			elif c2== "TG":
				inAr[4] += 1
			elif c2 in CH:
				inAr[3] += 1
			elif c2 in TH:
				inAr[7] += 1
			# 3 characters
			c3 = seq[i:i+3]
			if c3 in CHG:
				inAr[1] += 1
			elif c3 in THG:
				inAr[5] += 1
			elif c3 in CHH:
				inAr[2] += 1
			elif c3 in THH:
				inAr[6] += 1
		except IndexError:
			print('index error:', i)
	# end for
	return inAr

def primary( fastaFileStr ):
	cAr = readFasta(fastaFileStr)
	# print to output
	cArS = [ str(i) for i in cAr ]
	print( ' '.join(cArS) )
	
if __name__ == "__main__":
	if len(sys.argv) !=2:
		sys.stderr.write( "Usage: count_fasta.py <fast_file>" )
	else:
		fastaFileStr = sys.argv[1]
		primary(fastaFileStr)
