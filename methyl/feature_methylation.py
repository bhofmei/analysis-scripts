import sys, math, glob, multiprocessing, subprocess, os

# Usage: feature_methylation.py [-o=out_prefix] [-c=chrm_list] [-p=num_proc] [-m=meth_types] <-a | -g | -te> <gff_file> <allc_path> <sample> [sampleN]*

NUMPROC=1
CHRMLIST=["Chr1","Chr2","Chr3","Chr4","Chr5"]
MTYPES=["CG","CHG","CHH","C"]

def processInputs( gffFileStr, allcPath, sampleNamesAr, outPre, chrmList, numProc, mTypes, feature ):
	# feature: all, genes, te
	# generate dictionaries
	print( 'Reading GFF' )
	geneDict, teDict = readGFF( gffFileStr, feature )
	# process chromosomes
	print( 'Begin processing with {:d} processors'.format(numProc) )
	pool = multiprocessing.Pool(processes=numProc)
	results = [ pool.apply_async( processSample, args=( name, allcPath, chrmList, mTypes, geneDict, teDict ) ) for name in sampleNamesAr ]
	vals = [ p.get() for p in results ]
	# array of tuples (one tuple per sample) where tuple has array of weighted
	# methylation by context
	
	# write output
	outFileStr = outPre + "_" + feature + ".tsv"
	print( 'Writing output to', outFileStr )
	info = "#feature:{:s};chrms{:s};samples:{:s};meth:{:s}".format( feature, ','.join(chrmList), ','.join(sampleNamesAr), ','.join(mTypes) )
	writeOutput( outFileStr, vals, sampleNamesAr, mTypes, info )


def readGFF( gffFileStr, feature ):
	if feature == 'all':
		geneDict = {}
		teDict = {}
	elif feature == 'genes':
		geneDict = {}
		teDict = None
	elif feature == 'te':
		geneDict = None
		teDict = {}
	
	gffFile = open( gffFileStr, 'r' )
	for line in gffFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
		# (6) strand (7) ? (8) notes
		chrm = lineAr[0]
		start = int( lineAr[3] )
		end = int( lineAr[4] )
		# gene
		if lineAr[2] == 'gene' and feature in [ 'all', 'genes' ]:
			# check chrm
			if geneDict.get(chrm)==None:
				geneDict[chrm] = []
			geneDict[chrm] += [( start, end )]
		# te
		elif lineAr[2] in ['transposable_element', 'transposable_element_gene'] and feature in [ 'all', 'te' ]:
			if teDict.get(chrm)==None:
				teDict[chrm] = []
			teDict[chrm] += [( start, end )]
	# end for
	return geneDict, teDict

def processSample( name, allcPath, chrmList, mTypes, geneDict, teDict ):
	print( 'Processing', name )
	if geneDict != None:
		geneVals = [ [0,0] for x in mTypes ]
	else:
		geneVals = None
	if teDict != None:
		teVals = [ [0,0] for x in mTypes ]
	else:
		teVals = None
	
	# loop through chromosomes
	for chrm in chrmList:
		#print( teVals )
		tGene, tTE = processChrm( name, chrm, allcPath, mTypes, geneDict, teDict )
		if tGene != None:
			geneVals = addToArray( geneVals, tGene )
		if tTE != None:
			teVals = addToArray( teVals, tTE )
	# end for
	# compute wMeth -> array of float for mTypes
	wGene = computeMethylation( geneVals )
	wTE = computeMethylation ( teVals )
	return ( wGene, wTE )
		
def processChrm( name, chrm, allcPath, mTypes, geneDict, teDict ):
	
	#generate allc dictionary
	allcFileStr =os.path.normpath( '{:s}/allc_{:s}_{:s}.tsv'.format(allcPath, name, chrm ) )
	allcDict = readAllC( allcFileStr, mTypes )
	
	# compute for genes if necessary
	if geneDict == None:
		geneVals = None
	else:
		print( 'looking at genes for sample', name, 'chrm', chrm )
		geneVals = processFeature( allcDict, geneDict[chrm], mTypes )
		
	# compute for te if necessary
	if teDict == None:
		teVals = None
	else:
		print( 'looking at tes for sample', name, 'chrm', chrm )
		teVals = processChrmFeatures( allcDict, teDict[chrm], mTypes )
	
	# return results
	return geneVals, teVals

def readAllC( allCFileStr, mTypes ):
	allCFile = open( allCFileStr, 'r' )
	#print( 'Reading {:s}...'.format( allCFileStr ) )
	allCDict = {}
	
	for m in mTypes:
		allCDict[m] = {}
	
	for line in allCFile:
		if line.startswith( 'c' ):
			continue
		lineAr = line.rstrip().split('\t')
		# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
		# (6) methylated
		mLineType = findMethylType( lineAr[3] )
		if mLineType in mTypes or lineAr[3] in mTypes:
			allCDict[mLineType][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
		if 'C' in mTypes:
			allCDict['C'][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
	allCFile.close()
	return allCDict
	
def findMethylType( mc_class ):
	
	if mc_class.startswith( 'CG' ):
		return 'CG'
	elif mc_class.startswith( 'C' ) and mc_class.endswith( 'G' ):
		return 'CHG'
	else:
		return 'CHH'

def processChrmFeatures( allcDict, featureAr, mTypes ):
	
	mVal = [0] * len( mTypes )
	tVal = [0] * len( mTypes )
	
	# loop through mTypes
	for i in range(len(mTypes)):
		# loop through features
		for feat in featureAr:
			start = feat[0]
			end = feat[1]
			for pos in range( start, end+1 ):
				tup = allcDict.get(mTypes[i]).get( pos )
				if tup != None:
					mVal[i] += tup[0]
					tVal[i] += tup[1]
			# end for pos
		# end for feat
	# end for i
	z = zip( mVal, tVal )
	outAr = [(m,t) for m,t in z ]
	return outAr

def addToArray( oldArray, newArray ):
	# array has tuples/array of size 2 within
	
	for i in range(len(oldArray)):
		oldArray[i][0] += newArray[i][0]
		oldArray[i][1] += newArray[i][1]
	return oldArray

def computeMethylation( inAr ):
	if inAr == None:
		return None
	outAr = [-1] * len( inAr )
	for i in range(len(inAr)):
		if inAr[i][0] == 0 and inAr[i][1] == 0:
			outAr[i] = 0
		elif inAr[i][0] != 0 and inAr[i][1] == 0:
			print( methAr[i] )
		else:
			outAr[i] = float(inAr[i][0]) / float(inAr[i][1])
	return outAr
	
def writeOutput( outFileStr, valAr, sampleNamesAr, mTypes, info ):
	# valAr: sample -> tuple( feature ) -> array (mType) -> wMethylation
	header = "#sample\tfeature\tmType\twMethylation\n"
	outFile = open( outFileStr, 'w' )
	outFile.write( header + info + "\n" )
	
	# loop through samples
	for i in range(len(valAr)):
		# check features
		tup = valAr[i]
		if tup[0] != None:
			outFile.write( writeFeature( tup[0], sampleNamesAr[i], 'genes', mTypes ) )
		if tup[1] != None:
			outFile.write( writeFeature( tup[1], sampleNamesAr[i], 'tes', mTypes ) )
	# end for
	outFile.close()
			
def writeFeature( inAr, name, feature, mTypes ):
	outStr = ""
	# loop through array -> mTypes
	for i in range(len(mTypes)):
		outStr += '{:s}\t{:s}\t{:s}\t{:.6f}\n'.format( name, feature, mTypes[i], inAr[i] )
	return outStr

def parseInputs( argv ):
	outPre = 'out'
	chrmList = CHRMLIST
	mTypes = MTYPES
	numProc = NUMPROC
	startInd = 0
	feature = None
	
	for i in range(min(7,len(argv)-3)):
		if argv[i]=='-a':
			if feature != None:
				print( 'ERROR: do not specify more than one feature type' )
				exit()
			else:
				feature = 'all'
				startInd += 1
		elif argv[i] == '-g':
			if feature != None:
				print( 'ERROR: do not specify more than one feature type' )
				exit()
			else:
				feature = 'genes'
				startInd += 1
		elif argv[i] == '-te':
			if feature != None:
				print( 'ERROR: do not specify more than one feature type' )
				exit()
			else:
				feature = 'te'
				startInd += 1
		elif argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processers must be integer' )
				exit()
		elif argv[i].startswith( '-c=' ):
			chrmList = argv[i][3:].split(',')
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			mTmp = argv[i][3:].split(',')
			mTypes = [ x.upper() for x in mTmp ]
			startInd += 1
	# end for
	gffFileStr = argv[startInd]
	allcPath = argv[startInd+1]
	sampleNamesAr = []
	for j in range(startInd+2, len(argv)):
		sampleNamesAr += [ argv[j] ]
	
	processInputs( gffFileStr, allcPath, sampleNamesAr, outPre, chrmList, numProc, mTypes, feature )

if __name__ == "__main__":
	if len(sys.argv) < 5 :
		print ("Usage: feature_methylation.py [-o=out_prefix] [-c=chrm_list] [-p=num_proc] <-a | -g | -te> <gff_file> <allc_path> <sample> [sampleN]*")
	else:
		parseInputs( sys.argv[1:] )
