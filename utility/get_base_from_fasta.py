import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3.4 get_base_from_fasta.py <fasta_file>

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

def run( fastaDict ):
	print( "Enter 'exit' to quit " )
	print( "At the prompt, enter coordinate as\n  chromosome:position[-end_position]\nNote that positions are 1-based and chromosome must match chromosome titles from the fasta" )
	while True:
		txt = input( '> ' )
		if txt == 'exit':
			break
		ind1 = txt.find(':')
		if ind1 == -1:
			print( 'Input is not in the correct format' )
			continue
		chrm = txt[:ind1]
		fastaAr = fastaDict.get( chrm )
		if fastaAr == None:
			print( 'Chromosome specified does not exist for this FASTA' )
			continue
		pos = txt[ind1+1:]
		ind2 = pos.find('-')
		try:
			if ind2 == -1:
				print( '  ' + fastaAr[ int(pos) ] )
			else:
				start = int(pos[:ind2])
				end = int(pos[ind2+1:])
				print( '  ' + ''.join( fastaAr[start:end+1] ) )
		except ValueError:
			print( 'Position(s) are not proper integers' )
		except IndexError:
			print( 'Position(s) are out-of-bounds.  Length of {:s} is {:d}.'.format( chrm, len(fastaAr) ) )

if __name__ == "__main__":
	if len(sys.argv) != 2 :
		print ("Usage: python3.4 get_base_from_fasta.py <fasta_file>")
	else:
		fastaFileStr = sys.argv[1]
		print( 'Reading FASTA file...' )
		fastaDict = readFasta( fastaFileStr )
		run( fastaDict )
