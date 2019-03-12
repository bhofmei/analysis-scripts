import sys, math, glob, multiprocessing, subprocess, os
import numpy as np
import pandas as pd
import bth_util
from sklearn import linear_model
from functools import partial
from decodingpathcent import *
from centtransitions import Transitions

# Usage: python epigenotyping_pe_combbin_scaled.py [-u] [-c=bin_comb_thresh] [-d=decoding_type] [-p=num_proc] [-o=out_id] [-m=mother_sample] [-f=father_sample] [-b=bin_size] [-t=cent_start,cent_end] [-s=scale_factor] <input_file>

NUMPROC=1
BINSIZE=100000
DECODE='A'
UNIFORM=False
COMBINE=3
SCALE=1

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform, combineBins, cent, scaleFactor ):
	
	info = '#from_script: epigenotyping_pe_combbin_scaled.py; in_file:{:s}; bin_size:{:s}; decoding:{:s}; uni_class_prob:{:s}; combine_bins_threshold:{:d}; centromere:{:s}; scale_factor:{:g}'.format( os.path.basename( inFileStr ), bth_util.binSizeToStr( binSize ), formatDecoding( decoding).lower().replace('and',','), str(isUniform).lower(), combineBins, ('None' if cent == None else '{:s}-{:s}'.format( bth_util.binSizeToStr( cent[0] ), bth_util.binSizeToStr( cent[1] ) ) ), scaleFactor )
	print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
	print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	print( 'Uniform classification probabilities:', str(isUniform) )
	print( 'Decoding algorithm:', formatDecoding( decoding ) )
	print( 'Combine bin feature threshold:', combineBins )
	print( 'Centromere:', ( 'None' if cent == None else '{:s}-{:s}'.format( bth_util.binSizeToStr(cent[0]), bth_util.binSizeToStr(cent[1]) ) ) )
	print( 'Scale factor:', scaleFactor )
	
	# build dataframe
	print( ' Reading input file', os.path.basename( inFileStr ) )
	df = pd.read_table( inFileStr, header=1 )
	
	# check parent labels
	checkParents( df['sample'], parentLabelAr )
	
	# group by bin
	df['bin'] = df.pos // binSize
	transformation = None
	
	# get centromere bins if necessary
	if cent == None:
		centro = []
	else:
		cent = [ x // binSize for x in cent ]
		centro = list( range(cent[0], cent[1]+1) )
		#print(centro)
	
	# combine bins if necessary
	nbins = max(df['bin'])+1
	if combineBins > 0:
		print( ' Merging bins', end=' ... ' )
		df['tBin'] = df['bin']
		transformation = binTransformation( df, combineBins )
		# apply the transformation
		df['bin'] = df['tBin'].apply( lambda x: transformation[x] )
		
	
	dfBinGroup = df.groupby( 'bin' )
	if combineBins > 0:
		newNBins = len(dfBinGroup.groups )
		print( 'combined {:d} non-functional bins'.format( nbins - newNBins ) )
	
	# classify by bin
	print( ' Classifying {:d} bins with {:d} processors'.format( nbins, numProc ) )
	dfClass = runClassification( dfBinGroup, numProc, parentLabelAr, isUniform )
	dfClass.reset_index(inplace=True)
	#print( dfClass.head )
	del(df, dfBinGroup )
	# decode, if necessary
	if decoding != 'N':
		ignoreAr = parentLabelAr[:2] + ['MPV']
		transition = Transitions( dfClass, ignore = ignoreAr, cent=centro )
		transitionMatrix = transition.getTransitions()
		# write this matrix to file
		outFStr = determineTransFileName(inFileStr, outID, binSize, combineBins, scaleFactor )
		tLabels = [ 'mother', 'MPV', 'father' ]
		transData = pd.DataFrame( transitionMatrix, index=tLabels, columns= tLabels )
		with open( outFStr, 'w' ) as f:
			f.write(info+'\n')
		transData.to_csv( outFStr, sep='\t', mode='a' )
		
		# group by sample
		dfSampleGroup = dfClass.groupby( 'sample' )
		nsamples = len( dfSampleGroup.groups )
		
		print( ' {:s} decoding {:d} samples with {:d} processors'.format(  formatDecoding(decoding), nsamples, numProc ) )
		dfOutput = runDecoding( dfSampleGroup, numProc, transitionMatrix, decoding, centro, scaleFactor )
		dfOutput.set_index( ['bin', 'sample'], inplace=True )
		del( dfSampleGroup )
	else:
		dfOutput = dfClass
	
	# write output
	outFileStr = determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform, combineBins, scaleFactor )
	# if combination, undo transformation by applying the predictions to additional bins
	if combineBins > 0:
		dfOutput.reset_index(inplace=True)
		dfOutput['cBin'] = dfOutput['bin']
		dfOutputT = undoBinTransformation( dfOutput, transformation )
	else:
		dfOutputT = dfOutput.drop('cBin', axis=1)
	print( ' Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info+'\n')
	dfOutputT.to_csv( outFileStr, sep='\t', mode='a' )
	
	print( 'Done' )

def formatDecoding( decoding ):
	if decoding == 'V':
		return 'Viterbi'
	elif decoding == 'F':
		return 'Forward-backward'
	elif decoding == 'A':
		return 'Viterbi and Forward-backward'
	else:
		return 'None'

def checkParents( indexAr, parentLabelAr ):
	'''
		check parent labels
	'''
	#print(indexAr)
	#exit()
	# check mother
	i = np.where( indexAr == parentLabelAr[0] )
	if i[0].size == 0:
		print( 'ERROR: mother label {:s} not found'.format( parentLabelAr[0] ) )
		exit()
	# check father
	j = np.where( indexAr == parentLabelAr[1] )
	if j[0].size == 0:
		print( 'ERROR: father label {:s} not found'.format( parentLabelAr[1] ) )
		exit()

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
	

def runClassification( dfg, numProc, parentLabelAr, isUniform ):
	'''
		helper function that takes advantage of the number of processors
	'''
	mapfunc = partial( classifyBin, pla=parentLabelAr, u=isUniform )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat(res)

def classifyBin( group, pla=None, u=True ):
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

def classify( df, parentLabelAr, isUniform ):
	# rotate table
	dfs = df.pivot( index='sample', columns='pos', values='wei.meth' )
	# number of features used
	nfeat = dfs.shape[1]
	dfParent = dfs.loc[ [parentLabelAr[0], parentLabelAr[1]] ]
	# calculated mid-parent value
	mpv = dfParent.apply( np.mean, reduce=None )
	dfs.loc['MPV'] = mpv.transpose()
	
	# get training data
	train = dfs.loc[ [parentLabelAr[0], parentLabelAr[1], 'MPV'] ]
	trainClasses = np.array( train.index, dtype=np.str_ )
	trainClasses = renameParents( trainClasses, parentLabelAr )
	
	# create model
	if isUniform:
		clf = linear_model.LogisticRegression()
	else:
		clf = linear_model.LogisticRegression(class_weight={'mother':0.25,'MPV':0.5,'father':0.25})
	
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
	# replace mother
	if replacement == 1 or replacement == 3:
		i = np.where( inputAr == replaceAr[0] )
		inputAr[i[0][0]] = 'mother'
	if replacement == 2 or replacement == 3:
		i = np.where( inputAr == replaceAr[1] )
		inputAr[i[0][0]] = 'father'
	return inputAr

def runDecoding( dfg, numProc, transMat, decodeType, centro, scaleFactor ):
	mapfunc = partial( decodeSample, trans=transMat, d=decodeType, cent=centro, scale=scaleFactor )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat(res)

def decodeSample( group, trans=None, d='V', cent=[], scale=1 ):
	name, df = group
	#print(name)
	# get initial probabilities
	ip = getInitProb( df )
	res = decode( df, trans, ip, d, cent, scale )
	res['sample'] = name
	return res

def getInitProb( df ):
	labels = ['mother', 'MPV', 'father']
	counts = np.zeros( len(labels) )
	pred = df['prediction'].values
	for x in pred:
		y = labels.index(x)
		counts[y] += 1
	return counts / np.sum(counts)

def decode( df, transMat, initProb, decodeType, cent, scale ):
	# set up decoding
	if decodeType == 'A':
		alg = DecodeAll( df, transMat, initProb, cent, scale )
	elif decodeType == 'F':
		alg = DecodeForwardBackward( df, transMat, initProb, cent, scale )
	else:
		alg = DecodeViterbi( df, transMat, initProb, cent, scale )
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

def determineTransFileName( inFileStr, outID, binSize, combineBins, scaleFactor ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_epigenotype-trans_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_epigenotype-trans_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_epigenotype-trans_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
		
		# combining bins
	if combineBins > 0:
		outFileStr = outFileStr.replace('.tsv', '_cb-{:d}.tsv'.format( combineBins ) )
		
	if scaleFactor != 1:
		s = str( scaleFactor ).replace('.','-')
		outFileStr = outFileStr.replace( '.tsv', '_s{:s}.tsv'.format(s))
		
	return outFileStr

def determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform, combineBins, scaleFactor ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_epigenotype-results_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_epigenotype-results_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_epigenotype-results_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
		
	# combining bins
	if combineBins > 0:
		outFileStr = outFileStr.replace('.tsv', '_cb-{:d}.tsv'.format( combineBins ) )
	
	# scale factor
	if scaleFactor != 1:
		s = str( scaleFactor ).replace('.','-')
		outFileStr = outFileStr.replace( '.tsv', '_s{:s}.tsv'.format(s))
	
	# decoding and uniform
	if decoding != 'N' and isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni-{:s}.tsv'.format( 'vit' if decoding == 'V' else ('fb' if decoding == 'F' else 'vit-fb') ) )
	elif decoding != 'N':
		outFileStr = outFileStr.replace( '.tsv', '_{:s}.tsv'.format( 'vit' if decoding == 'V' else ('fb' if decoding == 'F' else 'vit-fb') ) )
	elif isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni.tsv' )
	
	return outFileStr
	
def parseInputs( argv ):
	numProc = NUMPROC
	binSize = BINSIZE
	outID = None
	parentLabelAr = ['mother', 'father', 0]
	decoding = DECODE
	isUniform = UNIFORM
	combineBins = COMBINE
	centromere=None
	scaleFactor = SCALE
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
				print( 'WARNING: number of processors must be integer...using default {:s}'.format(COMBINE) )
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
			else:
				print( 'WARNING: decoding option {:s} not recognized...using default {:s}'.format(opt, DECODE) )
			startInd += 1
		elif argv[i] == '-u':
			isUniform = True
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			tmp = argv[i][3:].split(',')
			tmp2 = [ bth_util.strToDistance( x ) for x in tmp ]
			if len(tmp2) != 2 or (False in tmp2):
				print( 'WARNING: centromere coordinates bad...not using' )
			else:
				centromere = tmp2
			startInd += 1
		elif argv[i].startswith( '-s=' ):
			try:
				scaleFactor = float( argv[i][3:] )
				startInd += 1
				if scaleFactor == 0:
					print( 'WARNING: scale factor must be greater than 0...using default', SCALE )
			except ValueError:
				print( 'WARNING: scale factor must be numeric...using default {:s}'.format(SCALE) )
				scaleFactor = SCALE		
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform, combineBins, centromere, scaleFactor )

def printHelp():
	print( 'Usage:\npython epigenotyping_pe_combbin_scaled.py [-u] [-c=bin_thresh] [-d=decoding_type][-p=num_proc]\n[-o=out_id] [-m=mother_sample] [-f=father_sample] [-b=bin_size] [-t=cent_start,cent_end] [-s=scale_factor] <input_file>' )
	print( 'Requried:' )
	print( 'input_file\tfile of of weighted methylation by position for samples' )
	print( 'Optional:' )
	print( '-u\t\tuniform class weights [default 1:2:1 for mother,\n\t\tMPV,father]' )
	print( '-d=decode_type\tdecoding type to use (capitlization ignored) [default {:s}]\n\t\tViterbi="v" or "viterbi"\n\t\tForward-Backward="forwardbackward", "f" or "fb"\n\t\tBoth="all" or "a"\n\t\tOff="false", "none", or "n"'.format(DECODE) )
	print( '-o=out_id\tidentifier for output file [default "out" or variation of\n\t\tinput file name]' )
	print( '-p=num_proc\tnumber of processors [default {:d}'.format(NUMPROC) )
	print( '-c=bin_thresh\tminimum number of features per bin to be classified\n\t\tgroups bins to reach this number [default {:d}'.format(COMBINE) )
	print( '-m=mother_label\tsample name of mother; for correct classification\n\t\t[default mother]' )
	print( '-f=father_label\tsample name of father; for correct classification\n\t\t[default father]' )
	print( '-b=bin_size\tsize of bins in bp [default {:s}]'.format( bth_util.binSizeToStr( BINSIZE ) ) )
	print( '-t=cent_start,cent_end\tcoordinates for centromere [default None]\n\t\twhen included, ignores this region for transitions and decoding' )
	print( '-s=scale_factor\tmultiplicative factor for weighting prediction probability over transition probability [default {:d} (unscaled)]'.format( SCALE ) )
	
if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
