import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import numpy as np
import pandas as pd
from sklearn import linear_model
import bth_util
from probpathsimple import ProbPathSimple
from simpletransitions import SimpleTransitions

# Usage: python epigenotype_by_logreg_pe.py [-n] [-u] [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_sample] [-f=father_sample] <input_file>

NUMPROC=1
BINSIZE=100000

def processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, isSmoothing, isUniform ):
	info = '#from_script:epigenotype_by_logreg.py; in_file:{:s}; bin_size:{:s}'.format( os.path.basename( inFileStr), bth_util.binSizeToStr( binSize ) )
	print( 'Weighted methylation file:', os.path.basename( inFileStr ) )
	print( 'Bin size:', bth_util.binSizeToStr( binSize ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	print( 'Smoothing:', str(isSmoothing) )
	print( 'Uniform classification probabilities:', str( isUniform ) )
	info += '; smoothing:{:s}; uni_class_prob:{:s}\n'.format( str(isSmoothing), str(isUniform) )

	# build data frame
	df = pd.read_table( inFileStr, header=1 )
	# check parent labels
	checkParents( df['sample'], parentLabelAr )

	# put in bins and analyze
	df['bin'] = df.pos // binSize
	nbins = max(df['bin'])+1
	dfg = df.groupby('bin')
	if numProc > 1:
		print( 'Begin classifying {:d} bins with {:d} processors'.format( nbins, numProc ) )
		res_class = runMultiprocessing( dfg, numProc, parentLabelAr, isUniform )
	else:
		print( 'Begin classifying {:d} bins'.format( nbins ) )
		res_class = dfg.apply( classLogRegImproved, pla=parentLabelAr, u=isUniform )
	res_class.reset_index(inplace=True)
	
	# smooth by sample
	if isSmoothing:
		ignoreAr = parentLabelAr + ['MPV']
		#transProbMat = computeTransitions( res_class, ignoreAr )
		transition = SimpleTransitions( res_class, ignore=ignoreAr )
		transProbMat = transition.run()
		print( transProbMat )
		groups = res_class.groupby( 'sample' )
		nsamples = len(groups.groups)
	
		# find optimum path for all samples, group by sample
		if numProc > 1:
			print( 'Begin smoothing {:d} samples with {:d} processors'.format(  nsamples, numProc ) )
			results = runMulti( groups, numProc, transProbMat )
		else:
			print( 'Begin smoothing {:d} samples'.format( nsamples ) )
			results = groups.apply( findOptimalPath, trans=transProbMat )
		results.set_index( ['bin', 'sample'], inplace=True )
	else:
		results = res_class
	
	# output file
	outFileStr = determineOutputFileName( inFileStr, outID, binSize, isSmoothing, isUniform )
	# write output
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info)
	results.to_csv( outFileStr, sep='\t', mode='a' )
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

def runMultiprocessing( dfg, numProc, parentLabelAr, isUniform ):
	from functools import partial
	mapfunc = partial( classLogRegImprovedMulti, pla=parentLabelAr, u=isUniform )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [group for group in dfg] )
	return pd.concat( res )

def classLogRegImproved( df, pla=None, u=False ):
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
	##### using log probabilities bc thats what we want for viterbi algorithm
	probs = clf.predict_log_proba( dfs )
	outdf = pd.DataFrame( probs, dfs.index, clf.classes_, None, True )
	outdf['prediction'] = pre
	return outdf

def classLogRegImprovedMulti( gr, pla=None, u=False ):
	name = gr[0]
	df = gr[1]
	rf = classLogRegImproved( df, pla, u )
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

def findOptimalPath( df, trans=None ):
	if trans == None:
		trans = [ [1.0/3.0]*3 for x in range(3) ]
	path = ProbPathSimple( df, trans )
	outDf = path.run()
	return outDf

def runMulti( groups, numProc, transitions ):
	from functools import partial
	mapfunc = partial( findOptimalPathMulti, trans=transitions )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [g for g in groups] )
	return pd.concat( res )

def findOptimalPathMulti( gr, trans=None ):
	name, df = gr
	res = findOptimalPath( df, trans )
	res['sample'] = name
	return res


def determineOutputFileName( inFileStr, outID, binSize, isSmoothing, isUniform ):
	outBaseName = os.path.basename( inFileStr )
	if outID == None:
		if '_wm_pos_' in inFileStr:
			outFileStr = outBaseName.replace( '_wm_pos_', '_logreg_{:s}_'.format( bth_util.binSizeToStr(binSize) ) )
		else:
			outFileStr = 'out_logreg_{:s}.tsv'.format( bth_util.binSizeToStr(binSize) )
	else:
		outFileStr = '{:s}_logreg_{:s}.tsv'.format( outID, bth_util.binSizeToStr(binSize) )
	if isSmoothing and isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni-opt.tsv' )
	elif isSmoothing:
		outFileStr = outFileStr.replace( '.tsv', '_opt.tsv' )
	elif isUniform:
		outFileStr = outFileStr.replace( '.tsv', '_uni.tsv' )	
	return outFileStr

def parseInputs( argv ):
	numProc = NUMPROC
	binSize = BINSIZE
	outID = None
	parentLabelAr = ['mother', 'father',0]
	isSmoothing = True
	isUniform = False
	startInd = 0

	for i in range(min(7,len(argv))):
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
		elif argv[i] == '-n':
			isSmoothing = False
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
	processInputs( inFileStr, numProc, binSize, outID, parentLabelAr, isSmoothing, isUniform )

def printHelp():
	print ("Usage: python epigenotype_by_logreg_pe.py [-n] [-p=num_proc] [-o=out_id] [-b=bin_size] [-m=mother_sample] [-f=father_sample] <input_file>")
	print( 'Requried:' )
	print( 'input_file\tfile of of weighted methylation by position for samples' )
	print( 'Optional:' )
	print( '-n\t\tturn off smoothing algorithm' )
	print( '-u\t\tuniform class weights [default 1:2:1 for mother,\n\t\tMPV,father]' )
	print( '-o=out_id\tidentifier for output file [default out or variation of\n\t\tinput file name]' )
	print( '-p=num_proc\tnumber of processors' )
	print( '-m=mother_label\tsample name of mother; for correct classification\n\t\t[default mother]' )
	print( '-f=father_label\tsample name of father; for correct classification\n\t\t[default father]' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
