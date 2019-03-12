import sys, glob, multiprocessing, subprocess, os, gzip

# Usage: python allc_to_home.py [-q] [-h] [-p=num_proc] [-o=out_dir] [-x=ex_chrms] <allc_file> [allc_file]*
# convert allc file to format needed by HOME DMR finding program

def processInputs( allcFileAr, outputDir, numProc, excludeList, isPrint ):
	if isPrint:
		print( 'Exclude chromosomes:', ', '.join(excludeList))
		print( 'Processing {:d} files with {:d} processors'.format(len(allcFileAr), numProc))
	
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(allcFile, outputDir, excludeList, isPrint) ) for allcFile in allcFileAr ]
	suc = [ p.get() for p in results ]

	print( 'Done' )

def processFile( allcFileStr, outputDir, excludeList, isPrint ):
	
	outFileStr = outFileName( allcFileStr, outputDir )
	outFile = open( outFileStr, 'w' )
	if allcFileStr.endswith('.gz'):
		allcFile = gzip.open(allcFileStr, 'rt' )
	else:
		allcFile = open( allcFileStr, 'r' )
			
	if isPrint:
		print('Converting {:s} to {:s}'.format(os.path.basename(allcFileStr), os.path.basename(outFileStr)))
	
	for line in allcFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		if line.startswith('#') or len(lineAr) < 7 or lineAr[6].isdigit() == False:
			continue
		if lineAr[0] in excludeList:
			continue
		mType = decodeMethType( lineAr[3] )
		# (0) chr1 (1) pos (2) strand (3) mType (4) mc_count (5) total
		if mType != '':
			outStr = '{:s}\t{:s}\t{:s}\n'.format('\t'.join(lineAr[0:3]), mType, '\t'.join(lineAr[4:6]))
			outFile.write(outStr)
	# end for line
	outFile.close()
	allcFile.close()
	return True

def decodeMethType( mStr ):

	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return ''
	else:
		return 'CHH'

def outFileName( allcFileStr, outputDir ):
	fileSplit = os.path.split(allcFileStr)
	outPath = (fileSplit[0] if outputDir == None else outputDir)
	
	l = fileSplit[1].replace('.gz', '' ).replace( '.gzip', '' )
	rInd = l.rfind( '.' )
	if rInd != -1:
		filePre = l[:rInd]
	else:
		filePre = l
	fileOut = filePre + '_home.tsv'
	return os.path.join(outPath, fileOut)
		
def parseInputs( argv ):
	isPrint = True
	numProc = 1
	outDir = None
	excludeList = []
	startInd = 0
	
	for i in range(min(5, len(argv)-1)):
		if argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i].startswith('-o='):
			outDir = argv[i][3:]
			startInd += 1
		elif argv[i].startswith('-x='):
			excludeList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be an integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for i
	
	allcFileAr = []
	for j in range(startInd, len(argv) ):
		allcFileAr += [ argv[j] ]
	
	processInputs( allcFileAr, outDir, numProc, excludeList, isPrint )

def printHelp():
	print( 'Usage:\tpython allc_to_home.py [-q] [-h] [-p=num_proc] [-o=out_dir] <allc_file> [allc_file]*')

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
