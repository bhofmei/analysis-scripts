import sys, math, glob, multiprocessing, subprocess, os
import pandas as pd
import numpy as np
import scipy.stats as stats
from functools import partial

# Usage: python compute_dmps_pe.py [-p=num_proc] [-o=out_id] [-v=min_cov] [-c=chrm_list] <allc_path> <sample1> <sample2> [sampleN]*
NUMPROC=1

def processInputs( allcPath, sampleNamesAr, outID, numProc, minCov, chrmList ):
	print( 'AllC path:', allcPath )
	print( 'Samples:', ', '.join(sampleNamesAr) )
	info = '#from_script: compute_dmps_pe.py'
	if minCov != None:
		info += '; min_cov: {:d}'.format( minCov )
	if chrmList != None:
		info += '; chrm_list: {:s}'.format(','.join(chrmList))
	info += '\n'
	
	# get allc files
	print( 'Reading {:d} allc files with {:d} processes'.format( len(sampleNamesAr), numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async(readAllc, args=(allcPath, sampleName, minCov) ) for sampleName in sampleNamesAr ]
	dfAr = [p.get() for p in results]
	
	# concat results
	print( 'Combining allc files' )
	dfCat = pd.concat(dfAr, ignore_index=True )
	del dfAr
	
	# if filter by chrm if necessary
	if chrmList != None:
		dfCat = dfCat[dfCat['chrm'].isin(chrmList)]
	
	# filter for methylated positions
	print( 'Filtering for positions methylated in 1+ samples' )
	dfCat.index = dfCat['sample']
	dfGroup = dfCat.groupby( ['chrm','pos'] )
	dfFilter = dfGroup.filter( lambda x: x['isMeth'].sum() > 0 )
	del dfCat, dfGroup
	
	# group by sample
	dfSample = dfFilter.groupby('sample')
	# get comparison list
	filterDict = getLabelFilter( list(dfSample.groups.keys()) )
	del dfFilter
	
	# call outer loop function
	print( 'Analyzing with {:d} processors'.format( numProc) )
	mapfunc = partial( runFisher, df=dfFilter, flt=filterDict )
	#dfFish = dfSample.apply( runFisher, df = dfFilter, flt = filterDict )
	res = pool.map( mapfunc, [g for g in dfSample] )
	dfFish = pd.concat( res )
	dfFish.set_index( 'sample_x', inplace=True )
	del dfSample
	
	# correct p-values
	print( 'Correcting p-values ({:d} tests)'.format( dfFish.shape[0] ) )
	dfFish['pvalue.adjust'] = pValueAdjust( dfFish['pvalue'] )
	
	# reorder columns
	dfOut = dfFish[ ['sample_y','chrm', 'pos', 'mCounts_x', 'tCounts_x', 'mCounts_y', 'tCounts_y', 'pvalue', 'pvalue.adjust']]
	
	del dfFish
	if outID == None:
		outID = 'out'
	outFileStr = outID + '_fisher_dmps.tsv'
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info)
	dfOut.to_csv( outFileStr, sep='\t', mode='a' )
	print( 'Done' )
	
def readAllc( allcPath, sampleName, minCov ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}{:s}.tsv'.format( allcPath, sampleName , ('' if minCov == None else '_cov'+str(minCov)) ) )
	print( ' {:s}'.format(os.path.basename(inFileStr)))
	df = pd.read_table( inFileStr,header=None,names=['chrm', 'pos', 'strand', 'mcClass','mCounts', 'tCounts', 'isMeth'],comment='#',usecols=['chrm','pos','mCounts','tCounts','isMeth'] )
	df['sample'] = sampleName
	return df

def getLabelFilter( labels ):
	filterDict = {}
	for i in range(len(labels)):
		filterDict[labels[i]] = []
		for j in range(i+1,len(labels)):
			filterDict[labels[i]] += [ labels[j] ]
		# end for j
	# end for i
	return filterDict

def runFisher( group, df, flt ):
	#curGroup = grDf.index.unique()[0]
	sample = group[0]
	grDf = group[1]
	fltList = flt[sample]
	dfMerge = pd.merge( grDf, df, on=['chrm','pos'])
	# get only comparisons that we want and drop unnecessary columns
	dfFilter = dfMerge.query( 'sample_y in @fltList' )
	dfFilterD = dfFilter.drop( ['sample_x','isMeth_x','isMeth_y'],axis=1 )

	# group and apply
	dfComp = dfFilterD.groupby( ['chrm', 'pos' ] )
	dfOutput = dfComp.apply( fisher )
	dfOutput['sample_x'] = sample
	return dfOutput

def fisher( data ):
	''' excutes fisher's exact test on incoming data
		returns original input data with additional column containing test
		p-value
	'''
	smp = data.loc[:,'sample_y'].iloc[0]
	mCounts = np.asarray(data.loc[:,['mCounts_x','mCounts_y']])[0]
	tCounts = np.asarray(data.loc[:,['tCounts_x','tCounts_y']])[0]
	subData = np.array([mCounts,tCounts])
	odr, p = stats.fisher_exact( subData )
	data['pvalue']=p
	return data
    
def pValueAdjust( pvalues ):
	'''	adjust p-values using BH
	'''
	n = len(pvalues)
	values = [ (p,i) for i,p in enumerate(pvalues) ]
	values.sort()
	values.reverse()
	newValues = []
	outValues = np.zeros(n)
	for i, vals in enumerate(values):
		rank = n - i
		pval,index = vals
		newValues.append( (n/rank)*pval )
	# end for i
	for i in range(n-1):
		if newValues[i] < newValues[i+1]:
			newValues[i+1] = newValues[i]
	# end for i
	for i, vals in enumerate(values):
		pval,index = vals
		outValues[index] = newValues[i]
	# end for i
	return outValues

def parseInputs( argv ):
	numProc = NUMPROC
	outID = None
	minCov=None
	chrmList = None
	startInd = 0

	for i in range(min(4,len(argv)-3)):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			tmpStr = argv[i][3:]
			chrmList = tmpStr.split(',')
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: number of processors must be integer...using 1' )
				numProc = NUMPROC
		elif argv[i].startswith( '-v=' ):
			try:
				minCov = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: minimum coverage must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	# end for
	allcPath = argv[startInd]
	if os.path.isdir( allcPath ) == False:
		print( 'ERROR: {:s} is not a path to a directory for allC files'.format( allcPath ) )
		exit()
	sampleNamesAr = []
	for i in range(startInd+1, len(argv)):
		sampleNamesAr += [ argv[i] ]
	processInputs( allcPath, sampleNamesAr, outID, numProc, minCov, chrmList )

def printHelp():
	print( 'Usage: python compute_dmps_pe.py [-p=num_proc] [-o=out_id] [-v=min_cov] [-c=chrm_list] <allc_path> <sample1> <sample2> [sampleN]*' )

if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
