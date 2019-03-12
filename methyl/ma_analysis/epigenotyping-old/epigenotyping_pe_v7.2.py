import sys, math, glob, multiprocessing, subprocess, os
import numpy as np
import pandas as pd
import bth_util
from sklearn import linear_model
from functools import partial
from decodingpath31 import *
from transitions import Transitions

# Usage: python epigenotyping_pe_v7 [-u | -e | -g=generation] [-q] [-c=bin_comb_thresh] [-d=decoding_type] [-p=num_proc] [-o=out_id] [-m=mother_samples] [-f=father_samples] [-b=bin_size] [-t=cent_start,cent_end[,cent2_start,cent2_end]] <input_file>

NUMPROC=1
BINSIZE=100000
DECODE='A'
GENERATION=2
ISPRINT=True
COMBINE=3

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, generation, combineBins, cent, isPrint ):
	
	info = '#from_script: epigenotyping_pe_v7.py; in_file:{:s}; bin_size:{:s}; decoding:{:s}; generation:{:d}; class_probs:{:s}; combine_bins_threshold:{:d}; centromere_{:s}'.format( os.path.basename( inFileStr ), bth_util.binSizeToStr( binSize ), formatDecoding( decoding).lower().replace('and',','), generation, ','.join( formatClassProbs( generation ) ), combineBins, ('None' if cent == None else '{:s}-{:s}'.format( bth_util.binSizeToStr( cent[0] ), bth_util.binSizeToStr( cent[1] ) ) ) )
	if isPrint:
		print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
		print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
		print( 'Mother label(s):', parentLabelAr[0] )
		print( 'Father label(s):', parentLabelAr[1] )
		print( 'Generation:', generation )
		print( 'Decoding algorithm:', formatDecoding( decoding ) )
		print( 'Combine bin feature threshold:', combineBins )
	if cent == None:
		centStr = 'None'
	else:
		centStr = ''
		for i in range(len(cent)//2):
			si = i*2
			centStr += '; {:s}-{:s}'.format( bth_util.binSizeToStr( cent[si] ), bth_util.binSizeToStr( cent[si+1] ) )
		centStr = centStr[2:]
			
	if isPrint:
		print( 'Centromere:', centStr )
	
	# build dataframe
	if isPrint:
		print( ' Reading input file', os.path.basename( inFileStr ) )
	df = pd.read_table( inFileStr, header=1 )
	
	# check parent labels
	newParentLabelAr = checkParents( df['sample'], parentLabelAr )
	tIgnoreAr = flattenList( newParentLabelAr[:2] )
	for i in range(len(newParentLabelAr[0])):
		tIgnoreAr += [ 'MPV{:d}'.format( i ) ]
	
	# group by bin
	df['bin'] = df.pos // binSize
	transformation = None
	
	# get centromere bins if necessary
	if cent == None:
		centBins = []
	else:
		cent = [ x // binSize for x in cent ]
		centBins = []
		#centBins = list( range(cent[0], cent[1]+1) )
		for i in range(len(cent) // 2 ):
			si = i * 2
			centBins += list( range(cent[si], cent[si+1]+1) )
	
	# combine bins if necessary
	nbins = max(df['bin'])+1
	if combineBins > 0:
		if isPrint:
			print( ' Merging bins', end=' ... ' )
		df['tBin'] = df['bin']
		transformation = binTransformation( df, combineBins )
		# apply the transformation
		df['bin'] = df['tBin'].apply( lambda x: transformation[x] )
		
	
	dfBinGroup = df.groupby( 'bin' )
	if combineBins > 0:
		newNBins = len( dfBinGroup.groups )
		info += '; non-functional_bins:{:d}'.format( nbins - newNBins )
		if isPrint:
			print( 'combined {:d} non-functional bins'.format( nbins - newNBins ) )
	
	# classify by bin
	if isPrint:
		print( ' Classifying {:d} bins with {:d} processors'.format( nbins, numProc ) )
	# get class probabilities based on generation number
	classProbs = computeClassProbs( generation )
	dfClass = runClassification( dfBinGroup, numProc, newParentLabelAr, classProbs )
	dfClass.reset_index(inplace=True)
	#print( dfClass.head )
	del(df, dfBinGroup )
	# decode, if necessary
	if decoding != 'N':
		transition = Transitions( dfClass, ignore = tIgnoreAr )
		transitionMatrix = transition.getTransitions()
		# write this matrix to file
		#outFStr = determineTransFileName(inFileStr, outID, binSize, combineBins )
		#tLabels = [ 'mother', 'MPV', 'father' ]
		#transData = pd.DataFrame( transitionMatrix, index=tLabels, columns= tLabels )
		#with open( outFStr, 'w' ) as f:
		#	f.write(info+'\n')
		#transData.to_csv( outFStr, sep='\t', mode='a' )
		
		# group by sample
		dfSampleGroup = dfClass.groupby( 'sample' )
		nsamples = len( dfSampleGroup.groups )
		tmpDecoding = ( 'F' if decoding == 'B' else decoding )
		if isPrint:
			print( ' {:s} decoding {:d} samples with {:d} processors'.format(  formatDecoding(tmpDecoding), nsamples, numProc ) )
		
		dfOutput = runDecoding( dfSampleGroup, numProc, transitionMatrix, tmpDecoding, centBins, classProbs )
		
		
		if decoding == 'B':
			dfNew = dfOutput.loc[:,['bin','sample']].copy()
			dfNew['MPV'] = np.log(dfOutput['fb.score.MPV'])
			dfNew['mother'] = np.log(dfOutput['fb.score.mother'])
			dfNew['father'] = np.log(dfOutput['fb.score.father'])
			dfNew['prediction'] = dfOutput['fb.prediction']

			transition = Transitions( dfNew, ignore = tIgnoreAr )
			transitionMatrix = transition.getTransitions()
			dfSampleGroup = dfNew.groupby( 'sample' )
			nsamples = len( dfSampleGroup.groups )
			
			if isPrint:
				print( ' {:s} decoding {:d} samples with {:d} processors'.format(  formatDecoding('V'), nsamples, numProc ) )
			dfOutputN = runDecoding( dfSampleGroup, numProc, transitionMatrix, 'V', centBins, classProbs )
			dfOutput[['vit.score.mother', 'vit.score.father', 'vit.score.MPV', 'vit.prob.mother', 'vit.prob.father', 'vit.prob.MPV', 'vit.prediction']] = dfOutputN[['vit.score.mother', 'vit.score.father', 'vit.score.MPV', 'vit.prob.mother', 'vit.prob.father', 'vit.prob.MPV', 'vit.prediction']]
			#print( dfOutput.head() )
		# end decoding == B
		dfOutput.set_index( ['bin', 'sample'], inplace=True )
		del( dfSampleGroup )
	else:
		dfOutput = dfClass
	
	# write output
	outFileStr = determineOutputFileName( inFileStr, outID, binSize, decoding, generation, combineBins )
	# if combination, undo transformation by applying the predictions to additional bins
	if combineBins > 0:
		dfOutput.reset_index(inplace=True)
		dfOutput['cBin'] = dfOutput['bin']
		dfOutputT = undoBinTransformation( dfOutput, transformation )
	else:
		dfOutputT = dfOutput.drop('cBin', axis=1)
	if isPrint:
		print( ' Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info+'\n')
	dfOutputT.to_csv( outFileStr, sep='\t', mode='a' )
	
	if isPrint:
		print( 'Done' )

def formatDecoding( decoding ):
	if decoding == 'V':
		return 'Viterbi'
	elif decoding == 'F':
		return 'Forward-backward'
	elif decoding == 'A':
		return 'Viterbi and Forward-backward'
	elif decoding == 'B':
		return 'Forward-backward then Viterbi'
	else:
		return 'None'

def computeClassProbs( n ):
	if n == 0:
		return None
	de = math.pow(2, n)
	nu1 = ( math.pow(2, n-1) ) - 1.0
	nu2 = 2
	oAr = [ x / de for x in [nu1, nu2, nu1] ]
	return oAr

def formatClassProbs (n):
	iAr = computeClassProbs( n )
	mAr = [ str(int(x * math.pow(2, n))) for x in iAr ]
	return mAr

def checkParents( indexAr, parentLabelAr ):
	'''
		check parent labels
	'''
	#print(indexAr)
	#exit()
	outMother = []
	outFather = []
	# check mother
	mList = parentLabelAr[0].split(',')
	for motherLabel in mList:
		i = np.where( indexAr == motherLabel )
		if i[0].size == 0:
			print( 'WARNING: mother label {:s} not found..excluding from analysis'.format( motherLabel ) )
		else:
			outMother += [ motherLabel ]
	# check father
	fList  = parentLabelAr[1].split(',')
	for fatherLabel in fList:
		j = np.where( indexAr == fatherLabel )
		if j[0].size == 0:
			print( 'WARNING: father label {:s} not found...excluding from analysis'.format( fatherLabel ) )
		else:
			outFather += [ fatherLabel ]
	
	# final checks
	if len( outMother ) == 0 or len( outFather ) == 0:
		print( 'ERROR: did not find samples for mother and/or father' )
		exit()
	elif len( outMother ) != len( outFather ):
		print( 'ERROR: number of mother and father samples not equal' )
		exit()
	else:
		return [outMother, outFather, parentLabelAr[2]]

def flattenList( inputList ):
	return [ item for sublist in inputList for item in sublist ]

def binTransformation( df, threshold ):
	# get dataframe for a single sample (it doesn't matter which)
	
	testSample = df.ix[0,'sample']
	#print( testSample )
	testData = df.loc[ df['sample'] == testSample,: ]
	testData = testData.copy()
	# count by bin
	counts = np.bincount( testData['bin'] )
	transformation = getTransformation( counts, threshold )
	return transformation
	# need to use the transformation to create a new column
	#return df['bin'].apply( lambda x: transformation[x] )

def getTransformation( countData, threshold ):
	nbins = len(countData)
	#print(nbins)
	transform = np.arange(nbins)
	mCounts = np.array(countData, copy=True)
	# loop backwards
	# count how many positions need change
	pCount = 0
	for i in range(nbins-1, 0, -1):
		# get count
		t = mCounts[i]
		# if there's enough, update transformation array
		if t >= threshold:
			for p in range(pCount):
				transform[i+p+1] -= p+1
			pCount = 0
		
		# if there's not enough, decrement transform[i] so it joins i-1 and
		# increase mCount of [i-1]
		else:
			pCount += 1
			mCounts[i-1] += mCounts[i]
	# end for
	return transform
	

def runClassification( dfg, numProc, parentLabelAr, classProbs ):
	'''
		helper function that takes advantage of the number of processors
	'''
	mapfunc = partial( classifyBin, pla=parentLabelAr, u=classProbs )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat(res)

def classifyBin( group, pla=None, u=[0.25, 0.5, 0.25] ):
	'''
		helper function that makes sure bin number is included in output
	'''
	name = group[0]
	df = group[1]
	res = classify( df, pla, u )
	res['bin'] = name
	res['sample']= res.index
	res.set_index( ['bin', 'sample'], inplace=True )
	return res

def classify( df, parentLabelAr, classProbs ):
	# rotate table
	dfs = df.pivot( index='sample', columns='pos', values='wei.meth' )
	# number of features used
	nfeat = dfs.shape[1]
	
	trainLabels = flattenList( parentLabelAr[:2] )
	nParents = len( parentLabelAr[0] )
	for i in range( nParents ):
		dfParent = dfs.loc[ [parentLabelAr[0][i], parentLabelAr[1][i]] ]
		# calculated mid-parent value
		mpv = dfParent.apply( np.mean, reduce=None )
		newLabel = 'MPV{:d}'.format( i )
		trainLabels += [ newLabel ]
		dfs.loc[newLabel] = mpv.transpose()
	
	# get training data
	train = dfs.loc[ trainLabels ]
	trainClasses = np.array( train.index, dtype=np.str_ )
	trainClasses = renameParents( trainClasses, parentLabelAr )
	
	# create model
	classWeights = {'mother': classProbs[0], 'MPV': classProbs[1], 'father': classProbs[2]}
	clf = linear_model.LogisticRegression(class_weight=classWeights)
	
	# fit training data
	clf.fit( train.values, trainClasses )
	# get classification probabilities and classification prediction
	prediction = clf.predict( dfs )
	logProbs = clf.predict_log_proba( dfs )
	# output data frame
	outdf = pd.DataFrame( logProbs, dfs.index, clf.classes_, None, True )
	outdf['prediction'] = prediction
	outdf['num.feat'] = nfeat
	del(dfs,train,clf)
	return outdf

def renameParents( inputAr, replaceAr ):
	replacement = replaceAr[2]
	nParents = len( replaceAr[0] )
	# replace mother
	if replacement == 1 or replacement == 3:
		for j in range( nParents ):
			i = np.where( inputAr == replaceAr[0][j] )
			inputAr[i[0][0]] = 'mother'
	# replace father
	if replacement == 2 or replacement == 3:
		for j in range( nParents ):
			i = np.where( inputAr == replaceAr[1][j] )
			inputAr[i[0][0]] = 'father'
	# replace MPV
	for j in range( nParents ):
		i = np.where( inputAr == 'MPV{:d}'.format(j) )
		inputAr[i[0][0]] = 'MPV'
	return inputAr

def runDecoding( dfg, numProc, transMat, decodeType, centro, classProbs ):
	mapfunc = partial( decodeSample, trans=transMat, d=decodeType, cent=centro, pi=classProbs )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat(res)

def decodeSample( group, trans=None, d='V', cent=[], pi=None ):
	name, df = group
	#print( name )
	res = decode( df, trans, d, cent, pi )
	res['sample'] = name
	return res

def decode( df, transMat, decodeType, cent, pi ):
	# set up decoding
	if decodeType == 'A':
		alg = DecodeAll( df, transMat, cent, pi )
	elif decodeType == 'F':
		alg = DecodeForwardBackward( df, transMat, cent, pi )
	else:
		alg = DecodeViterbi( df, transMat, cent, pi )
	# run decoding
	outdf = alg.run()
	return outdf
	
def undoBinTransformation( df, transformation ):
	# get new dictionary
	transDict = undoTransformation( transformation )
	dfCopy = df.copy()
	#dfCopy.reset_index(inplace=True)
	
	# loop through keys in trans dict
	for key in transDict.keys():
		# get all rows with cBin == key
		dfNew = df[ df['cBin'] == key ]
		# loop through appending new rows as necessary
		binReplacements = transDict[key]
		for index,row in dfNew.iterrows():
			for b in binReplacements:
				row.loc['bin'] = b
				dfCopy = dfCopy.append(row, ignore_index=True)
		# end for index
	# end for key
	# reset index and return
	dfCopy.sort_values( by=['sample','bin'], inplace=True )
	dfCopy.set_index( ['bin', 'sample'], inplace=True )
	
	return dfCopy

def undoTransformation( transformation ):
	outDict = {}
	for i in range(len(transformation)):
		t = transformation[i]
		if i != t: # bin was replaced
			# if not in dict, add
			if outDict.get( t ) == None:
				outDict[t] = []
			outDict[t] += [ i ]
	# end for i
	return outDict

def determineTransFileName( inFileStr, outID, binSize, combineBins ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_transition_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_transition_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_transition_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
		
		# combining bins
	if combineBins > 0:
		outFileStr = outFileStr.replace('.tsv', '_cb-{:d}.tsv'.format( combineBins ) )
		
		return outFileStr

def determineOutputFileName( inFileStr, outID, binSize, decoding, generation, combineBins ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_epigenotype-v7.2_{:s}_g-{:d}_'.format( bth_util.binSizeToStr(binSize), generation ) )
		else:
			outFileStr = 'out_epigenotype-v7.2_{:s}_g-{:d}.tsv'.format( bth_util.binSizeToStr(binSize), generation )
	else:
		outFileStr = '{:s}_epigenotype-v7.2_{:s}_g-{:d}.tsv'.format( outID, bth_util.binSizeToStr(binSize), generation )
		
	# combining bins
	if combineBins > 0:
		outFileStr = outFileStr.replace('.tsv', '_cb-{:d}.tsv'.format( combineBins ) )
	
	# decoding
	if decoding != 'N':
		outFileStr = outFileStr.replace( '.tsv', '_{:s}.tsv'.format( 'vit' if decoding == 'V' else ('fb' if decoding == 'F' else ('both' if decoding == 'B' else 'vit-fb')) ) )

	return outFileStr
	
def parseInputs( argv ):
	numProc = NUMPROC
	binSize = BINSIZE
	outID = None
	parentLabelAr = ['mother', 'father', 0]
	decoding = DECODE
	generation = GENERATION
	combineBins = COMBINE
	isPrint = ISPRINT
	centromere=None
	startInd = 0
	
	for i in range( min(9, len(argv)-1) ):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-b=' ):
			inStr = argv[i][3:]
			binSize = bth_util.strToDistance( inStr )
			if binSize == False:
				print( 'WARNING: cannot convert {:s} to bin size...using default {:s}'.format( inStr, bth_util.binSizeToStr(BINSIZE) ) )
				binSize = BINSIZE
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be integer...using 1' )
				numProc = NUMPROC
		elif argv[i].startswith( '-c=' ):
			try:
				combineBins = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: combine bins must be integer...using default {:s}'.format(COMBINE) )
				combineBins = COMBINE
		elif argv[i].startswith( '-m=' ):
			parentLabelAr[0] = argv[i][3:]
			parentLabelAr[2] += 1
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			parentLabelAr[1] = argv[i][3:]
			parentLabelAr[2] += 2
			startInd += 1
		elif argv[i].startswith( '-d=' ):
			opt = argv[i][3:].lower()
			if opt == 'false' or opt == 'none' or opt== 'n':
				decoding = 'N'
			elif opt == 'viterbi' or opt == 'v':
				decoding = 'V';
			elif opt == 'forwardbackward' or opt == 'f' or opt == 'fb':
				decoding = 'F'
			elif opt == 'all' or opt == 'a':
				decoding = 'A'
			elif opt == 'both' or opt == 'b':
				decoding = 'B'
			else:
				print( 'WARNING: decoding option {:s} not recognized...using default viterbi'.format(opt) )
			startInd += 1
		elif argv[i].startswith( '-g=' ):
			try:
				generation = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: generation must be integer...using default {:s}'.format(COMBINE) )
				generation = COMBINE
		elif argv[i] == '-q':
			isPrint = False
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			tmp = argv[i][3:].split(',')
			tmp2 = [ bth_util.strToDistance( x ) for x in tmp ]
			if len(tmp2) % 2 != 0 or (False in tmp2):
				print( 'WARNING: centromere coordinates bad...not using' )
			else:
				centromere = tmp2
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, generation, combineBins, centromere, isPrint )

def printHelp():
	print( 'Usage:\npython epigenotyping_pe_v7.2.py [-q] [-g=generation]  [-c=bin_thresh] [-d=decoding_type][-p=num_proc]\n[-o=out_id] [-m=mother_sample(s)] [-f=father_sample(s)] [-b=bin_size] <input_file>' )
	print( 'Requried:' )
	print( 'input_file\tfile of of weighted methylation by position for samples' )
	print( 'Optional:' )
	print( '-q\t\tquiet, do not print progress' )
	print( '-g=generation\tgeneration of self-crossing; used to determine\n\t\tclassification probabilities; use 0 for uniform weight\n\t\t[default 2]' )
	print( '-d=decode_type\tdecoding type to use (capitlization ignored) [default {:s}]\n\t\tViterbi="v" or "viterbi"\n\t\tForward-Backward="forwardbackward", "f" or "fb"\n\t\tBoth="all" or "a"\n\t\tOff="false", "none", or "n"'.format(DECODE) )
	print( '-o=out_id\tidentifier for output file [default "out" or variation of\n\t\tinput file name]' )
	print( '-p=num_proc\tnumber of processors [default {:d}'.format(NUMPROC) )
	print( '-c=bin_thresh\tminimum number of features per bin to be classified\n\t\tgroups bins to reach this number [default {:d}'.format(COMBINE) )
	print( '-m=mother_samples\tsample name(s) of mother; for correct classification\n\t\t[default mother]' )
	print( '-f=father_samples\tsample name(s) of father; for correct classification\n\t\t[default father]' )
	print( '-b=bin_size\tsize of bins in bp [default {:s}]'.format( bth_util.binSizeToStr( BINSIZE ) ) )
	
if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
