import sys, math, glob, multiprocessing, subprocess, os

# Usage: allC_to_annoj.py [-o=out_dir] [-c=chrm_list] [-p=num_proc] <allC_path> <sample_name> [sample_name]*

NUMPROC=2
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

# assembly | position | strand | class | mc | h 
# only include methylated positions

def processInputs( allCPath, sampleNamesAr, outDir, chrmList, numProc ):
	
	# check outDir
	if outDir == None:
		outDir = allCPath
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processSample, args=(sample, allCPath, chrmList, outDir) ) for sample in sampleNamesAr ]
	nums = [ p.get() for p in results ]
	print( 'Done.' )

def processSample( sampleName, allCPath, chrmList, outDir ):
	
	for chrm in chrmList:
		allCFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, sampleName, chrm ) )
		outFileStr = os.path.normpath( '{:s}/mc_{:s}_{:s}.tsv'.format( outDir, sampleName.replace('-','_'), chrm.replace('Chr','') ) )
		print( 'Processing {:s} {:s}...'.format( sampleName, chrm ) )
		readAllC( allCFileStr, outFileStr )
	return 1

def readAllC( allCFileStr, outFileStr ):
	allCFile = open( allCFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if int( lineAr[6] ):
			lineAr[0] = lineAr[0].replace('Chr','').replace('chr','')
			lineAr[3] = decodeMethType( lineAr[3] )
			if lineAr[3] == False:
				continue
			outStr = '{:s}\n'.format( '\t'.join( lineAr[0:6] ) )
			outFile.write( outStr )
	# end for
	allCFile.close()
	outFile.close()

def parseInputs( argv ):
	# allC_to_annoj.py [-o=out_dir] [-c=chrm_list] [-p=num_proc] <allC_path> <sample_name> [sample_name]*
	
	chrmList = CHRMLIST
	numProc = NUMPROC
	outDir = None
	startInd = 0
	
	for i in range( min(3,len(argv)) ):
		if argv[i].startswith( '-o=' ):
			outDir = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i].startswith( '-c=' ):
			chrmList = argv[i][3:].split(',')
			startInd += 1
	# end for
	
	allCPath = argv[startInd]
	if os.path.isdir( allCPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory'.format(allCPath) )
			exit()
			
	sampleNamesAr = []
	for j in range( startInd+1, len(argv) ):
		sampleNamesAr += [ argv[j] ]
	
	processInputs( allCPath, sampleNamesAr, outDir, chrmList, numProc )
				

def decodeMethType( mStr ):
	
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return False
	else:
		return 'CHH'
	
		

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: allC_to_annoj.py [-o=out_dir] [-c=chrm_list] [-p=num_proc] <allC_path> <sample_name> [sample_name]*")
	else:
		parseInputs( sys.argv[1:] )
