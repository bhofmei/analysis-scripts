import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import numpy as np
import pandas as pd
from sklearn import linear_model
import bth_util

# Usage: python epigenotype_by_logreg_pe.py [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_sample] [-f=father_sample] <input_file>

NUMPROC=1
BINSIZE=100000

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr ):
	info = '#from_script:epigenotype_by_logreg.py; in_file:{:s}; bin_size:{:s}'.format( os.path.basename( inFileStr), bth_util.binSizeToStr( binSize ) )
	print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
	print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )

	# build data frame
	df = pd.read_table( inFileStr, header=1 )
	# check parent labels
	checkParents( df['sample'], parentLabelAr )

	# put in bins and analyze
	df['bin'] = df.pos // binSize
	nbins = max(df['bin'])+1
	dfg = df.groupby('bin')
	if numProc > 1:
		print( 'Begin processing {:d} bins with {:d} processors'.format( nbins, numProc ) )
		results = runMultiprocessing( dfg, numProc, parentLabelAr )
	else:
		print( 'Begin processing {:d} bins'.format( nbins ) )
		results = dfg.apply( classLogRegImproved, pla=parentLabelAr )
	print( 'Finished processing' )
	# output file
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_logreg_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_logreg_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_logreg_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
	# write output
	print( 'Writing output to', outFileStr )
	results.to_csv( outFileStr, sep='\t' )
	print( 'Done' )

def checkParents( indexAr, parentLabelAr ):
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

def runMultiprocessing( dfg, numProc, parentLabelAr ):
	from functools import partial
	mapfunc = partial( classLogRegImprovedMulti, pla=parentLabelAr )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat( res )

def classLogRegImproved( df, pla=None ):
	# get MDV
	dfs = df.pivot(index='sample',columns='pos',values='wei.meth')
	dfsm = dfs.loc[[pla[0], pla[1]]]
	mdv = dfsm.apply(np.mean,reduce=None)
	dfs.loc['MDV'] = mdv.transpose()
	clf = linear_model.LogisticRegression()
	# get training set
	train = dfs.loc[[pla[0], pla[1], 'MDV']]
	tr_cl = np.array(train.index,dtype=np.str_)
	tr_cl = renameParents( tr_cl, pla )
	# create model and fit initial data
	clf = linear_model.LogisticRegression()
	clf.fit( train.values, tr_cl )
	# get the class probabilities
	pre = clf.predict( dfs )
	#probs = clf.predict_proba( dfs )
	##### using log probabilities bc thats what we want for viterbi algorithm
	probs = clf.predict_log_proba( dfs )
	#dec = clf.decision_function( dfs )
	outdf = pd.DataFrame( probs, dfs.index, clf.classes_, None, True )
	outdf['prediction'] = pre
	#outdf['sample'] = dfs.index
	return outdf

def classLogRegImprovedMulti( gr, pla=None ):
	name = gr[0]
	df = gr[1]
	rf = classLogRegImproved( df, pla )
	rf['bin'] = name
	rf['sample'] = rf.index
	rf.set_index( ['bin','sample'],inplace=True )
	return rf

def renameParents( inAr, replaceAr ):
	replacement = replaceAr[2]
	# replace mother
	if replacement == 1 or replacement == 3:
		i = np.where( inAr == replaceAr[0] )
		inAr[i[0][0]] = 'mother'
	if replacement == 2 or replacement == 3:
		i = np.where( inAr == replaceAr[1] )
		inAr[i[0][0]] = 'father'
	return inAr

def parseInputs( argv ):
	numProc = NUMPROC
	binSize = BINSIZE
	outID = None
	parentLabelAr = ['mother', 'father',0]
	startInd = 0

	for i in range(min(5,len(argv))):
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
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr )

def printHelp():
	print ("Usage: python epigenotype_by_logreg_pe.py [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_sample] [-f=father_sample] <input_file>")


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
