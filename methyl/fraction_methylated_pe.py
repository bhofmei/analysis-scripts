import sys, math, glob, multiprocessing, subprocess, os
'''
	reports the number of methylated, unmethylated, and low coverage positions
	for input samples along with weighted methylation of all positions and
	subsets
'''
# Usage: python3.4 fraction_methylated_pe.py [-o=out_prefix] [-m=meth_types] [-c=chromosomes] [-v=min_cov] [-p=num_proc] <allC_path> <sample_name> [sample_name]*
NUMPROC = 1
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def processInputs( allCPath, lineNamesAr, outPre, minCov, chrmList, mTypes, numProc ):

	outFileStr = outPre + '_cov{:d}.tsv'.format( minCov )
	
	print( 'Methylation types: {:s}\nChromosomes: {:s}\nMinimum coverage: {:d}\nLines included: {:s}\nOutput written to: {:s}\n'.format( ' '.join(mTypes), ' '.join( chrmList ), minCov, ' '.join(lineNamesAr), outFileStr ) )
	
	totalDict = {}
	for m in mTypes:
		totalDict[m] = [ [0]*9 for x in lineNamesAr ]
	print( 'Begin processing with {:d} processors...'.format( numProc ) )
	# multiprocess
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processChrm, args=(x, allCPath, lineNamesAr, minCov, mTypes)) for x in chrmList ]
	aDicts = [ p.get() for p in results ]
	
	# combine the analyzed dicts
	print( 'Combining analyzed dictionaries...' )
	for d in aDicts:
		totalDict = combineAnalyzeDicts( totalDict, d )
	
	# output
	info = '#from_script:fraction_methylated_pe.py; min_coverage:{:d}; chrom_included:{:s}'.format( minCov, ','.join( chrmList ) )
	print( 'Writing output...' ) 
	writeOutput( outFileStr, totalDict, lineNamesAr, info )
	print( 'Done.' )

def processChrm( chrm, allCPath, lineNamesAr, minCov, mTypes ):
	
	print( 'Starting on {:s}...'.format( chrm ) )
	# generate chrmDict
	chrmDict = generateChrmDict( chrm, allCPath, lineNamesAr, minCov, mTypes )
	# analyze chrmDict
	print( 'Analyzing {:s}...'.format( chrm ) )
	analyzeDict = analyzeChrmDict( chrmDict, lineNamesAr )
	print( 'Finished with {:s}.'.format( chrm ) )
	return analyzeDict

def generateChrmDict( chrm, allCPath, lineNamesAr, minCov, mTypes ):
	
	chrmDict = {}
	for m in mTypes:
		chrmDict[m] = {}
	
	# loop through lines
	for i in range(len(lineNamesAr)):
		allCFileStr = os.path.normpath('{:s}/allc_{:s}_{:s}.tsv'.format( allCPath, lineNamesAr[i], chrm ))
		# get allC diction
		print( 'Reading {:s}...'.format( allCFileStr ) )
		tmpDict = readAllC( allCFileStr, minCov, mTypes )
		# add to chrmDict
		chrmDict = addToChrmDict( chrmDict, tmpDict, i )
	return chrmDict

def readAllC( allCFileStr, minCov, mTypes ):
	
	allCFile = open( allCFileStr, 'r' )
	
	allCDict = {}
	for m in mTypes:
		allCDict[m] = {}
	for line in allCFile:
		# header
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		m = determineMethylType( lineAr[3] )
		pos = int( lineAr[1] )
		methCount = int( lineAr[4] )
		totalCount = int( lineAr[5] )
		methylated = int( lineAr[6] )
		tup = (methCount, totalCount, (methylated if totalCount >= minCov else -1 ) )
		if m in mTypes:
			allCDict[m][pos] = tup
		if 'C' in mTypes:
			allCDict['C'][pos] = tup
	# end for line
	allCFile.close()
	return allCDict

def determineMethylType( mStr ):
	
	if mStr.startswith( 'CG' ):
		return 'CG'
	elif mStr.startswith('C') and mStr.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'
		
def addToChrmDict( chrmDict, tmpDict, lineNum ):
	
	mTypes = chrmDict.keys()
	for m in mTypes:
		chrmDict[m] = combineDicts( chrmDict[m], tmpDict[m], lineNum )
	return chrmDict
	
def combineDicts( chrmDict, tmpDict, lineNum ):
	
	# loop through all positions in chrmDict
	for pos in chrmDict.keys():
		tup = tmpDict.get( pos )
		# not in tmpDict -> add (0,0,-1) tuple
		if tup == None:
			chrmDict[pos] += [ (0,0,-1) ]
		# in tmpDict -> add to chrmDict, remove pos from tmpDict
		else:
			chrmDict[pos] += [ tup ]
			del tmpDict[pos]
	# loop through remaining positions in tmpDict
	for pos in tmpDict.keys():
		prevAr = [ (0,0,-1) for x in range(lineNum) ]
		chrmDict[pos] = prevAr + [ tmpDict[pos] ]
	
	return chrmDict

def analyzeChrmDict( chrmDict, lineNamesAr ):
	analyzeDict = {}
	
	for m in chrmDict.keys():	# mTypes
		#print(m)
		analyzeDict[m] = analyzeChrmDict_M( chrmDict[m], lineNamesAr )
	return analyzeDict
	
def analyzeChrmDict_M( chrmDict, lineNamesAr ):
	# stats: (0) num_pos_meth (1) num_pos_unmeth (2) num_pos_na 
	# (3) num_mc-read_meth (4) num_tot-read_meth (5) num_mc-read_meth 
	# (6) num_tot-read_meth (7) num_mc-read (8) num_tot-red
	outMatr = [ [0]*9 for x in lineNamesAr ]
	
	# loop through positions
	for pos in chrmDict.keys():
		# loop through lines
		tupAr = chrmDict[pos]
		for i in range(len( lineNamesAr ) ):
			# methylated
			if tupAr[i][2] == 1:
				outMatr[i][0] += 1
				# add methylated read counts
				outMatr[i][3] += tupAr[i][0]
				outMatr[i][4] += tupAr[i][1]
				# add to total read counts
				outMatr[i][7] += tupAr[i][0]
				outMatr[i][8] += tupAr[i][1]
			# unmethylated
			elif tupAr[i][2] == 0:
				outMatr[i][1] += 1
				# add methylated read counts
				outMatr[i][5] += tupAr[i][0]
				outMatr[i][6] += tupAr[i][1]
				# add to total read counts
				outMatr[i][7] += tupAr[i][0]
				outMatr[i][8] += tupAr[i][1]
			# unknown
			elif tupAr[i][2] == -1:
				#print( lineNamesAr[i]  + '\t' + str(pos) )
				outMatr[i][2] += 1
		# end for lines
	# end for positions
	return outMatr
	
def combineAnalyzeDicts( totalDict, tmpDict ):
	
	#loop through mtypes
	for m in totalDict.keys():
		# loop through lines
		mtDict = totalDict[m]
		mDict = tmpDict[m]
		# loop through lines
		for j in range(len( mtDict ) ):
			mtDict[j] = [ mtDict[j][i] + mDict[j][i] for i in range(len(mtDict[j])) ]
	return totalDict

def writeOutput( outFileStr, totalDict, lineNamesAr, info ):
	
	outFile = open( outFileStr, 'w' )
	# header - mType line num_pos_meth num_pos_unmeth num_pos_na weighted_meth_meth weighted_meth_unmeth weighted_meth_overall
	header = '{:s}\nmeth.type\tline.name\tnum.pos.meth\tnum.pos.unmeth\tnum.pos.na\tweighted.meth\tweighted.unmeth\tweighted.overall\n'.format( info )
	outFile.write( header )
	
	# loop through meth types
	for methType in totalDict.keys():
		outMatrix = totalDict[methType]
		# loop through lines
		for j in range(len(outMatrix)):
			intStrAr = [ '{:d}'.format( outMatrix[j][i]) for i in range(3) ]
			floatStrAr = [ '{:.4f}'.format( float(outMatrix[j][i]) / float(outMatrix[j][i+1]) * 100)for i in range(3,8,2) ]
			outStr = methType + '\t' + lineNamesAr[j] + '\t' + '\t'.join( intStrAr ) + '\t' + '\t'.join( floatStrAr ) + '\n'
			outFile.write( outStr )
		# end for lines
	# end for meth types
	outFile.close()

def parseInputs( argv ):
	
	mTypes = ['C','CG','CHG','CHH']
	outPre = 'frac_pos_meth'
	minCov = 0
	numProc = NUMPROC
	chrmList = CHRMLIST
	startInd = 0
	
	for i in range(min(len(argv)-2, 5) ):
		if argv[i].startswith( '-o='):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m='):
			mTypes = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-c='):
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-v='):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum coverage must be integer')
				exit()
		elif argv[i].startswith( '-p='):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer')
				exit()
	# end for
	
	allCPath = argv[startInd]
	if os.path.isdir( allCPath ) == False:
			print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allCPath ) )
			exit()
	
	lineNamesAr = []
	for i in range( startInd + 1, len(argv) ):
		lineNamesAr += [ argv[i] ]
	
	processInputs( allCPath, lineNamesAr, outPre, minCov, chrmList, mTypes, numProc)
	
if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python3.4 fraction_methylated_pe.py [-o=out_prefix] [-m=meth_types] [-c=chromosomes] [-v=min_cov] [-p=num_proc] <allC_path> <sample_name> [sample_name]*")
	else:
		parseInputs( sys.argv[1:] )
