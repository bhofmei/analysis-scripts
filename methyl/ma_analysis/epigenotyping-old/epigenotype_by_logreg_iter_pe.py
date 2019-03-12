import sys, math, glob, multiprocessing, subprocess, os
import numpy as np
import pandas as pd
from sklearn import linear_model
import bth_util
from transitions import Transitions
from probpathsimple import ProbPathSimple
from probpath import *


# Usage: python epigenotype_by_logreg_iter_pe.py [-u] [-s=decode_type] [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_sample] [-f=father_sample] [-x=max_iter] <input_file>

NUMPROC = 1
BINSIZE=100000
DECODE='V'
UNIFORM=False
MAXITER=10

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform, maxIter ):
	dType = ('Viterbi' if decoding == 'V' else ('Forward-backward' if decoding == 'F' else ('Viterbi and Forward-backward' if decoding == 'A' else 'None') ) )
	info = '#from_script:epigenotype_by_logreg.py; in_file:{:s}; bin_size:{:s}; decoding:{:s}; uni_class_prob:{:s}\n'.format( os.path.basename( inFileStr), bth_util.binSizeToStr( binSize ), dType.lower().replace(' and ', ','), str(isUniform) )
	print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
	print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	print( 'Uniform classification probabilities:', str( isUniform ) )
	print( 'Decoding algorithm:', dType)
	
	# build data frame
	df = pd.read_table( inFileStr, header=1 )
	# check parent labels
	checkParents( df['sample'], parentLabelAr )
	
	# group by bin and analyze
	df['bin'] = df.pos // binSize
	nbins = max(df['bin'])+1
	dfg = df.groupby('bin')
	if numProc > 1:
		print( 'Begin classifying {:d} bins with {:d} processors'.format( nbins, numProc ) )
		res_class = runMultiClassification( dfg, numProc, parentLabelAr, isUniform )
	else:
		print( 'Begin classifying {:d} bins'.format( nbins ) )
		res_class = dfg.apply( classLogReg, pla=parentLabelAr, u=isUniform )
	res_class.reset_index(inplace=True)
	
	# decode if necessary
	if decoding != 'N':
		ignoreAr = parentLabelAr + ['MPV']
		print( 'Generating transition matrix' )
		transition = Transitions( res_class, ignore=ignoreAr )
		transitions = transition.getTransitions()
		# find optimum path for all samples
		groups = res_class.groupby( 'sample' )
		nsamples = len(groups.groups)
		
		if numProc > 1:
			print( 'Begin {:s} decoding {:d} samples with {:d} processors'.format(  dType, nsamples, numProc ) )
			results = runMultiPath( groups, numProc, transitions, isUniform, decoding )
		else:
			print( 'Begin {:s} decoding {:d} samples'.format( dType, nsamples ) )
			results = groups.apply( findOptimalPath, trans=transitions, u=isUniform, d=decoding )
		results.set_index( ['bin', 'sample'], inplace=True )
	else:
		results = res_class
	
	# output file
	outFileStr = determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform )
	# write output
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info)
	results.to_csv( outFileStr, sep='\t', mode='a' )
	print( 'Done' )

def checkParents( indexAr, parentLabelAr ):
	'''
		checks that the parent labels are part of the samples
		if not, exits the program
	'''
	#indexAr = np.array(indexObj,dtype=np.str_)
	# check mother
	i = np.where( indexAr == parentLabelAr[0] )
	if i[0].size == 0:
		print( 'ERROR: mother label {:s} not found'.format( parentLabelAr[0] ) )
		print( 'Samples:', indexAr )
		exit()
	j = np.where( indexAr == parentLabelAr[1] )
	if j[0].size == 0:
		print( 'ERROR: father label {:s} not found'.format( parentLabelAr[1] ) )
		print( 'Samples:', indexAr )
		exit()

def classLogReg( df, pla=None, u=False ):
	'''
		df is input data frame for all samples in one bin
		pla is parent label array
		u indicates if uniform class weights should be use
		l indicates if probabilities should be log probs
	'''
	# get MDV
	dfs = df.pivot(index='sample',columns='pos',values='wei.meth')
	dfsm = dfs.loc[[pla[0], pla[1]]]
	mdv = dfsm.apply(np.mean,reduce=None)
	dfs.loc['MPV'] = mdv.transpose()
	clf = linear_model.LogisticRegression()
	
	# get training set
	train = dfs.loc[[pla[0], pla[1], 'MPV']]
	tr_cl = np.array(train.index,dtype=np.str_)
	tr_cl = renameParents( tr_cl, pla )
	# create model and fit initial data
	if u:
		clf = linear_model.LogisticRegression()
	else:
		clf = linear_model.LogisticRegression(class_weight={'mother':0.25,'MPV':0.5,'father':0.25})
	clf.fit( train.values, tr_cl )
	# get the class probabilities
	pre = clf.predict( dfs )
	# use log probabilities
	probs = clf.predict_log_proba( dfs )
	outdf = pd.DataFrame( probs, dfs.index, clf.classes_, None, True )
	outdf['prediction'] = pre
	return outdf	
	
def runMultiClassification( dfg, numProc, parentLabelAr, isUniform ):
	'''
		main helper function for running classification with
		multiple processors
		necessary because classification algorithm has additional parameters
	'''
	from functools import partial
	mapfunc = partial( classLogRegMulti, pla=parentLabelAr, u=isUniform )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat( res )

def classLogRegMulti( gr, pla=None, u=False ):
	'''
		helper for when using more than one processor
		makes sure resulting data frame is the same as if run with single 
		processor
		gr is a group--tuple of label and data frame
	'''
	name = gr[0]
	df = gr[1]
	rf = classLogReg( df, pla, u )
	rf['bin'] = name
	rf['sample'] = rf.index
	rf.set_index( ['bin','sample'],inplace=True )
	return rf

def renameParents( inAr, replaceAr ):
	'''
		rename parent classifiers based on replaceAr
	'''
	replacement = replaceAr[2]
	# replace mother
	if replacement == 1 or replacement == 3:
		i = np.where( inAr == replaceAr[0] )
		inAr[i[0][0]] = 'mother'
	if replacement == 2 or replacement == 3:
		i = np.where( inAr == replaceAr[1] )
		inAr[i[0][0]] = 'father'
	return inAr

def findOptimalPath( df, trans=np.array([]), u=False, d='V' ):
	if trans.size == 0:
		trans = np.full( (3,3),1.0/3 )
	if d == 'A':
		path = ProbPathAll( df, trans, u )
	elif d == 'F':
		path = ProbPathFB( df, trans, u )
	else:
		#path = ProbPathViterbi( df, trans, u )
		path = ProbPathSimple( df, trans)
	outDf = path.run()
	return outDf

def runMultiPath( groups, numProc, transitions, uniform, decoding ):
	from functools import partial
	mapfunc = partial( findOptimalPathMulti, trans=transitions, u=uniform, d=decoding )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [g for g in groups] )
	return pd.concat( res )

def findOptimalPathMulti( gr, trans=None, u=False, d='V' ):
	name, df = gr
	res = findOptimalPath( df, trans, u, d )
	res['sample'] = name
	return res

def getCrossovers	

def determineOutputFileName( inFileStr, outID, binSize, decoding, isUniform ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_logreg_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_logreg_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_logreg_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
	
	if decoding != 'N' and isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni-{:s}.tsv'.format( 'vit' if decoding == 'V' else ('fb' if decoding == 'F' else 'vit_fb') ) )
	elif decoding != 'N':
		outFileStr = outFileStr.replace( '.tsv', '_{:s}.tsv'.format( 'vit' if decoding == 'V' else ('fb' if decoding == 'F' else 'vit_fb') ) )
	elif isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni.tsv' )	
	return outFileStr
	
def parseInputs( argv ):
	numProc = NUMPROC
	binSize = BINSIZE
	outID = None
	parentLabelAr = ['mother', 'father',0]
	decoding = DECODE
	isUniform = UNIFORM
	maxIter = MAXITER
	startInd = 0

	for i in range(min(8,len(argv))):
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
		elif argv[i].startswith( '-x=' ):
			try:
				maxIter = int( argv[i][3:] )
				startInt += 1
			except ValueError:
				print( 'WARNING: maximum iterations must be integer...using 10' )
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, decoding, isUniform, maxIter )

def printHelp():
	print( 'Usage: python epigenotype_by_logreg_pe.py [-u] [-d=decode_type] [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_samples] [-f=father_samples] <input_file>' )
	print( 'Requried:' )
	print( 'input_file\tfile of of weighted methylation by position for samples' )
	print( 'Optional:' )
	print( '-u\t\tuniform class weights [default 1:2:1 for mother,\n\t\tMPV,father]' )
	print( '-d=decode_type\tdecoding type to use (capitlization ignored) [default V]\n\t\tViterbi="v" or "viterbi"\n\t\tForward-Backward="forwardbackward", "f" or "fb"\n\t\tBoth="all" or "a"\n\t\tOff="false", "none", or "n"' )
	print( '-o=out_id\tidentifier for output file [default out or variation of\n\t\tinput file name]' )
	print( '-p=num_proc\tnumber of processors' )
	print( '-m=mother_label\tsample name of mother; for correct classification\n\t\t[default mother]' )
	print( '-f=father_label\tsample name of father; for correct classification\n\t\t[default father]' )
	print( '-x=max_iter\tmaximum iterations for determining crossovers [default 10]' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
