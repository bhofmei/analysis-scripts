import sys, multiprocessing, subprocess, os, gzip

# Usage: python find_mpos_diff_pe.py [-v=min_cov] [-c=chrm_list] [-o=out_id] [-p=num_proc] [-m=meth_types] <allc_path> <sample1_name> <sample2_name> 
# get individual positions that differ based on binomial test
NUMPROC = 1
MINCOV=3
CHRMLIST = [ 'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5' ]

def processInputs( allcPath, sample1name, sample2name, numProc, outID, methType, minCov, chrmList, isPrint ):
	if isPrint:
		print( 'AllC path:', allcPath )
		print( 'Samples:', sample1name, '&', sample2name )
		print( 'Methylation types:',  ', '.join(methType) )
		print( 'Minimum coverage:', minCov )
	print( 'Chromosomes:', ', '.join(chrmList) )
	info = '#from_script: find_mpos_diff_pe.py; samples: {:s}; methyl_types: {:s}; min_cov: {:d}; chrms: {:s}'.format( ','.join( [sample1name, sample2name] ), ','.join(methType), minCov, ','.join(chrmList) )
	
	# initial check for files
	suc1 = checkFiles( allcPath, sample1name, chrmList )
	suc2 = checkFiles( allcPath, sample2name, chrmList )
	if not(suc1 and suc2):
		exit()
	
	# process by chrm
	if isPrint:
		print( 'Begin processing by chrm with', numProc, 'processors' )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(processChrm, args=(allcPath, x, sample1name, sample2name, minCov, methType, isPrint) ) for x in chrmList ]
	outArs = [p.get() for p in results]
	
	outFileStr = '{:s}_mpos_diff.tsv'.format( outID )
	if isPrint:
		print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, outArs, info )
	print( 'Done' )
	

def checkFiles( allcPath, sample, chrmList ):
	
	for chrm in chrmList:
		if os.path.isfile( os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sample, chrm ) ) ) == False and os.path.isfile( os.path.normpath('{:s}/allc_{:s}_{:s}.tsv.gz'.format( allcPath, sample, chrm ) ) ) == False:
			print( 'ERROR: allC file for sample {:s} chrm {:s} not found'.format( sample, chrm ) )
			return False
	return True
	
def processChrm( allcPath, chrm, sample1name, sample2name, minCov, methType, isPrint ):
	if isPrint:
		print( 'Processing', chrm )
	outDict = {0: chrm}
	
	allcFile1 = os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sample1name, chrm ) )
	isZip1 = False
	if os.path.isfile( allcFile1 ) == False:
		allcFile1 = allcFile1 + '.gz'
		isZip1 = True
	
	allcFile2 = os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format( allcPath, sample2name, chrm ) )
	isZip2 = False
	if os.path.isfile( allcFile2 ) == False:
		allcFile2 = allcFile2 + '.gz'
		isZip2 = True
	
	# read file 1
	outDict = readAllc( allcFile1, minCov, methType, outDict, isZip1 )
	# read file 2
	outDict = readAllc( allcFile2, minCov, methType, outDict, isZip2 )
	# filter dictionary
	outStr = filterDict( outDict )
	del outDict
	
	return outStr

def readAllc( allcFileStr, minCov, methTypes, inDict, isZip ):
	
	if isZip:
		allCFile = gzip.open( allcFileStr, 'rt' )
	else:
		allCFile = open( allcFileStr, 'r' )
	for line in allCFile:
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		# headers and min coverage check
		if len(lineAr) < 7 or lineAr[6].isdigit() == False or int( lineAr[5] ) < minCov:
			continue
		mType = decodeMethType( lineAr[3] )
		if mType not in methTypes:
			continue
		pos = int( lineAr[1] )
		check = inDict.get( pos )
		if check == None:
			inDict[pos] = []
		inDict[pos] += [ int(lineAr[6]) ]
	# end for
	allCFile.close()
	return inDict

def decodeMethType( mStr ):
	
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.endswith( 'G' ):
		return 'CHG'
	elif mStr == 'CNN':
		return 'CNN'
	else:
		return 'CHH'

def filterDict( inDict ):
	chrm = inDict.pop( 0, 'NA' )
	outStr = ''
	for pos in sorted( inDict.keys() ):
		a = inDict[pos]
		# only keep when len = 2 and sum = 1
		if len( a ) == 2 and sum( a ) == 1:
			outStr += '{:s}\t{:d}\n'.format( chrm, pos )
			#outStr += '{:s}\t{:d}\t{:d}\t{:d}\n'.format( chrm, pos, a[0], a[1] )
	# end for
	return outStr

def writeOutput( outFileStr, outAr, info ):
	outFile = open( outFileStr, 'w' )
	outFile.write( info + '\nchrm\tpos\n' )
	for s in sorted(outAr):
		outFile.write( s )
	outFile.close()
	

def parseInputs( argv ):
	numProc = 1
	outID = 'out'
	methType = ['CG','CHG','CHH']
	minCov = MINCOV
	chrmList = CHRMLIST
	isPrint = True
	startInd = 0
	
	for i in range(min(6,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			m = argv[i][3:].split(',')
			tmpM = []
			for x in m:
				if x.upper() not in ['CG','CHG','CHH']:
					print( 'WARNING: {:s} is invalid methylation type..ignoring' )
				else:
					tmpM += [x.upper()]
			methType = ( methType if len(tmpM) == 0 else tmpM )
			startInd += 1
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be integer..using default', NUMPROC )
				startInd += 1
				numProc = NUMPROC
		elif argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: minimum coverage must be integer...using default', MINCOV )
				minCov = MINCOV
				startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	allcPath = argv[startInd]
	sample1name = argv[startInd + 1]
	sample2name = argv[startInd + 2]
	processInputs( allcPath, sample1name, sample2name, numProc, outID, methType, minCov, chrmList, isPrint )
	

def printHelp():
	print( 'Usage:\tpython find_mpos_diff_pe.py [-h] [-q] [-v=min_cov] [-c=chrm_list]\n\t[-o=out_id] [-p=num_proc] [-m=meth_types] <allc_path> <sample1_name>\n\t<sample2_name>' )
	print()
	print( 'Required:' )
	print( 'allc_path\tpath to allc files' )
	print( 'sample_name\tnames of samples to compare' )
	print()
	print( 'Optional:' )
	print( '-h\t\tprint this help screen and exit' )
	print( '-q\t\tquiet; does not print progress' )
	print( '-v=min_cov\tmin coverage to include a position [default {:d}]'.format( MINCOV ) )
	print( '-o=out_id\tstring for output file name [default "out"]' )
	print( '-c=chrm_list\tcomma-separated list of chrms [default arabidopsis]' )
	print( '-p=num_proc\tnum processors to use [default {:d}]'.format( NUMPROC ) )
	print( '-m=meth_types\tcomma-separated list of "CG", "CHG", and/or "CHH"\n\t\t[default all]' )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
