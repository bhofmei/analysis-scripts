import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import numpy as np
import pandas as pd
from probpath import ProbPath

# Usage: python epigenotype_smoothing_pe.py [-p=numProc] [-m=mother_label] [-f=father_label] <input_file>

NUMPROC=1

def processInputs( inFileStr, numProc, parentLabelAr ):
	
	print( 'Input file:', os.path.basename( inFileStr ) )
	print( 'Mother label:', parentLabelAr[0] )
	print( 'Father label:', parentLabelAr[1] )
	# read in file
	df = pd.read_table( inFileStr )
	# check parent labels
	checkParents( df['sample'], parentLabelAr )
	
	# get transition probabilities
	# ignore MDV and parents
	ignoreAr = parentLabelAr + ['MDV']
	transProbMat = computeTransitions( df, ignoreAr )
	_printTransMat( transProbMat )
	# using class probabilities and transition probabilities, 
	# find optimum path for all samples
	# group by sample
	groups = df.groupby( 'sample' )
	nsamples = len(groups.groups)
	
	if numProc > 1:
		print( 'Begin processing {:d} samples with {:d} processors'.format(  nsamples, numProc ) )
		results = runMultiprocessing( groups, numProc, transProbMat )
	else:
		print( 'Begin processing {:d} samples'.format( nsamples ) )
		results = groups.apply( findOptimalPath, trans=transProbMat )
	print( 'Finished processing' )
	results.set_index( ['bin', 'sample'], inplace=True )
	# output results
	rInd = inFileStr.rfind( '.' )
	outFileStr = inFileStr[:rInd] + '_opt' + inFileStr[rInd:]
	print( 'Writing output to', outFileStr )
	results.to_csv( outFileStr, sep='\t' )
	print( 'Done' )

def checkParents( indexAr, parentLabelAr ):
	#indexAr = np.array(indexObj,dtype=np.str_)
	# check mother
	i = np.where( indexAr == parentLabelAr[0] )
	if i[0].size == 0:
		print( 'WARNING: mother label {:s} not found'.format( parentLabelAr[0] ) )
	j = np.where( indexAr == parentLabelAr[1] )
	if j[0].size == 0:
		print( 'WARNING: father label {:s} not found'.format( parentLabelAr[1] ) )

def computeTransitions( df, ignoreLabelsAr ):
	# nine transition types
	labels = ['mother', 'MDV','father','total']
	transCountDict = {}
	for x in labels:
		transCountDict[x] = {}
		for y in labels:
			transCountDict[x][y] = 1
	dfp = df.pivot( index='sample', columns = 'bin', values='prediction' )
	for index,row in dfp.iterrows():
		if index not in ignoreLabelsAr:
			for i in range(row.size - 1):
				try:
					tmpDict = transCountDict[row.iloc[i]]
					tmpDict[row.iloc[i]] += 1
					tmpDict['total'] += 1
				except KeyError:
					print( 'key', row.iloc[i], 'or', row.iloc[i+1], 'not found' )
	# end for index, row
	transProbMat = []
	for x in labels:
		ar = [ math.log(float(transCountDict[x][y]) / transCountDict[x]['total']) for y in labels ]
		transProbMat.append( ar )
	return transProbMat

def findOptimalPath( df, trans=None ):
	if trans == None:
		trans = [ [1.0/3.0]*3 for x in range(3) ]
	path = ProbPath( df, trans )
	outDf = path.run()
	return outDf

def runMulti( groups, numProc, transitions ):
	from functools import partial
	mapfunc = partial( fileOptimalPathMulti, trans=transitions )
	with multiprocessing.Pool(processes=numProc) as p:
		res = p.map( mapfunc, [g for g in groups] )

def findOptimalPathMulti( gr, trans=None ):
	name, df = gr
	res = findOptimalPath( df, trans )
	res['bin']  = res.index
	res['sample'] = name
	return res

def _printTransMat( inMat ):
	outStr = ''
	lAr = ['m','h','f']
	headerAr = [ ' {:>6s}'.format( x ) for x in lAr ]
	outStr += ' '*2 + ''.join(headerAr) + '\n' + ' '*3 + '-'*(3*7) +'\n'
	for i in range(len(lAr)):
		tmpStr = ' {:s}|'.format( lAr[i] )
		transAr = [ ' {:>6.4f}'.format(inMat[i][j]) for j in range(3) ]
		outStr += tmpStr + ''.join(transAr)+ '|\n'
	outStr += ' '*3 + '-'*(3*7)
	print( outStr )
	
		
def parseInputs( argv ):
	numProc = NUMPROC
	parentLabelAr = ['mother','father']
	startInd = 0
	
	for i in range(min(3,len(argv))):
		if argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print('WARNING: number of processors must be integer...using 1')
				numProc = NUMPROC
		elif argv[i].startswith( '-m=' ):
			parentLabelAr[0] = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			parentLabelAr[1] = argv[i][3:]
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	inFileStr = argv[startInd]
	processInputs( inFileStr, numProc, parentLabelAr )
	

def printHelp():
	print(  'Usage: python epigenotype_smoothing_pe.py [-p=numProc] [-m=mother_label] [-f=father_label] <input_file>' )
	print( 'Requried:' )
	print( 'input_file\tfile of computed probabilities for bins with the predictions\n\t\tshould be the output of epigenotype by logreg' )
	print( 'Optional:' )
	print( '-p=num_proc\tnumber of processors to use for decoding optimal path' )
	print( '-m=mother_label\tsample name of mother; to exclude from transition calculations' )
	print( '-f=father_label\tsample name of father' )
	

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
