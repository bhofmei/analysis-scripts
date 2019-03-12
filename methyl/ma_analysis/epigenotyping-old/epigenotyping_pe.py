import sys, math, glob, multiprocessing, subprocess, os
import numpy as np
import pandas as pd
import bth_util
from sklearn import linear_model
from functools import partial
from decodingpath import *
from transitions import Transitions

# Usage: python epigenotyping_pe.py [-u] [-d=decoding_type] [-p=num_proc] [-o=out_id] [-m=mother_sample] [-f=father_sample] [-b=bin_size] <input_file>

NUMPROC=1
BINSIZE=100000
DECODE='A'
UNIFORM=False

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform ):
	
	info = '#from_script: epigenotyping_pe.py; in_file:{:s}; bin_size:{:s}; decoding:{:s}; uni_class_prob:{:s}'.format( os.path.basename( inFileStr ), bth_util.binSizeToStr( binSize ), formatDecoding( decoding).lower().replace('and',','), str(isUniform).lower() )
	print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
	print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	print( 'Uniform classification probabilities:', str(isUniform) )
	print( 'Decoding algorithm:', formatDecoding( decoding ) )
	
	# build dataframe
	df = pd.read_table( inFileStr, header=1 )
	
	# check parent labels
	checkParents( df['sample'], parentLabelAr )
	
	# group by bin
	df['bin'] = df.pos // binSize
	nbins = max(df['bin'])+1
	dfBinGroup = df.groupby( 'bin' )
	
	# classify by bin
	print( 'Begin classifying {:d} bins with {:d} processors'.format( nbins, numProc ) )
	dfClass = runClassification( dfBinGroup, numProc, parentLabelAr, isUniform )
	dfClass.reset_index(inplace=True)
	#print( dfClass.head )
	del(df, dfBinGroup )
	# decode, if necessary
	if decoding != 'N':
		ignoreAr = parentLabelAr[:2] + ['MPV']
		transition = Transitions( dfClass, ignore = ignoreAr )
		transitionMatrix = transition.getTransitions()
		
		# group by sample
		dfSampleGroup = dfClass.groupby( 'sample' )
		nsamples = len(dfSampleGroup.groups )
		
		print( 'Begin {:s} decoding {:d} samples with {:d} processors'.format(  formatDecoding(decoding), nsamples, numProc ) )
		dfOutput = runDecoding( dfSampleGroup, numProc, transitionMatrix, decoding )
		dfOutput.set_index( ['bin', 'sample'], inplace=True )
		del( dfSampleGroup )
	else:
		dfOutput = dfClass
	
	# write output
	outFileStr = determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform )
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info+'\n')
	dfOutput.to_csv( outFileStr, sep='\t', mode='a' )
	
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

def runDecoding( dfg, numProc, transMat, decodeType ):
	mapfunc = partial( decodeSample, trans=transMat, d=decodeType )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat(res)

def decodeSample( group, trans=None, d='V' ):
	name, df = group
	res = decode( df, trans, d )
	res['sample'] = name
	return res

def decode( df, transMat, decodeType ):
	# set up decoding
	if decodeType == 'A':
		alg = DecodeAll( df, transMat )
	elif decodeType == 'F':
		alg = DecodeForwardBackward( df, transMat )
	else:
		alg = DecodeViterbi( df, transMat )
	# run decoding
	outdf = alg.run()
	return outdf

def determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_epigenotype_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_epigenotype_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_epigenotype_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
	
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
	startInd = 0
	
	for i in range( min(6, len(argv)-1) ):
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
				print( 'WARNING: decoding option {:s} not recognized...using default viterbi'.format(opt) )
			startInd += 1
		elif argv[i] == '-u':
			isUniform = True
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform )

def printHelp():
	print( 'Usage: python epigenotyping_pe.py [-u] [-d=decoding_type] [-p=num_proc] [-o=out_id] [-m=mother_sample] [-f=father_sample] [-b=bin_size] <input_file>' )
	print( 'Requried:' )
	print( 'input_file\tfile of of weighted methylation by position for samples' )
	print( 'Optional:' )
	print( '-u\t\tuniform class weights [default 1:2:1 for mother,\n\t\tMPV,father]' )
	print( '-d=decode_type\tdecoding type to use (capitlization ignored) [default A]\n\t\tViterbi="v" or "viterbi"\n\t\tForward-Backward="forwardbackward", "f" or "fb"\n\t\tBoth="all" or "a"\n\t\tOff="false", "none", or "n"' )
	print( '-o=out_id\tidentifier for output file [default out or variation of\n\t\tinput file name]' )
	print( '-p=num_proc\tnumber of processors' )
	print( '-m=mother_label\tsample name of mother; for correct classification\n\t\t[default mother]' )
	print( '-f=father_label\tsample name of father; for correct classification\n\t\t[default father]' )
	print( '-b=bin_size\tsize of bins in bp [default {:s}]'.format( bth_util.binSizeToStr( BINSIZE ) ) )
	
if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
