import sys, math, glob, multiprocessing, subprocess, os, bisect, random
def warn( *args, **kwargs ):
	pass
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.warn = warn
import pandas as pd
import numpy as np
from sklearn import metrics

# Usage: simulation_accuracy.py [-q] [-o=out_id] <input_file>

def processInputs( inFileStr, outID, isQuiet ):

	if not isQuiet:
		print( 'Reading input file' )
	# get information
	tbl = pd.read_csv( inFileStr )
	
	# get class specific scores
	if not isQuiet:
		print( 'Getting global scores' )
	otbl = getMetrics( tbl )
	otbl['sample'] = 'global'
	otbl['test'] = 'global'
	
	if not isQuiet:
		print( 'Getting sample scores' )
	grouptbl = tbl.groupby('sample')
	gtbl = grouptbl.apply(getMetrics)
	# remove unnecessary level
	gtbl.index = gtbl.index.droplevel(level=1)
	gtbl.reset_index(inplace=True)
	gtbl['test'] = 'global'
	
	if not isQuiet:
		print( 'Getting test scores' )
	testtbl = tbl.groupby('test')
	ttbl = testtbl.apply(getMetrics)
	ttbl.index = ttbl.index.droplevel(level=1)
	#ttbl.index.names = ['sample']
	ttbl.reset_index(inplace=True)
	ttbl['sample'] = 'global'
	
	if not isQuiet:
		print( 'Getting individual scores' )
	indivtbl = tbl.groupby(['test', 'sample'])
	itbl = indivtbl.apply(getMetrics)
	# remove unnecessary level
	itbl.index = itbl.index.droplevel(level=2)
	itbl.reset_index(inplace=True)
	
	if not isQuiet:
		print( 'Combining' )
	outtbl = pd.concat([otbl, gtbl, ttbl, itbl], ignore_index=True)
	
	if outID == None:
		bn = os.path.basename( inFileStr )
		rInd = bn.rfind('.')
		if rInd != -1:
			outID = bn[:rInd]
		else:
			outID = bn
	outFileStr = outID + '_metrics.csv'
	
	with open( outFileStr, 'w' ) as f:
		f.write('#from_script: simulation_accuracy.py; in_file: {:s}\n'.format(os.path.basename(inFileStr)) )
	outtbl.to_csv( outFileStr, mode='a' )
	
	
def getMetrics( tbl ):
	lbl = ['m', 'h', 'f']
	cPrec, cRecall, cF1, cSup = metrics.precision_recall_fscore_support(tbl.expected, tbl.prediction, labels=lbl, warn_for=())
	#cAcc = np.array([-1]*3)
	gPrec, gRecall, gF1, gSup = metrics.precision_recall_fscore_support(tbl.expected, tbl.prediction, labels=lbl, average='micro',warn_for=())
	#gAcc = metrics.jaccard_similarity_score(tbl.expected, tbl.prediction)
	prec = np.append(cPrec, gPrec)
	rec = np.append(cRecall, gRecall)
	f1 = np.append(cF1, gF1)
	#acc = np.append( cAcc, gAcc )
	labels = lbl + ['global']
	
	#outtbl = pd.DataFrame( { 'label': labels, 'precision': prec, 'recall': rec, 'f1-score': f1, 'accuracy': acc } )
	outtbl = pd.DataFrame( { 'label': labels, 'precision': prec, 'recall': rec, 'f1-score': f1 } )
	return outtbl

def parseInputs( argv ):
	outID = None
	isQuiet = False
	startInd = 0
	
	for i in range(min(2, len(argv)-1)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i] == '-q':
			isQuiet = True
			startInd += 1
	# end for
	inFileStr = argv[startInd]
	processInputs( inFileStr, outID, isQuiet )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python simulation_accuracy.py [-q] [-o=out_id] <input_file>")
	else:
		parseInputs( sys.argv[1:] )
