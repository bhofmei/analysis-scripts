import sys, math, glob, multiprocessing, subprocess, os
import pandas as pd
import numpy as np
import scipy.stats as stats
from functools import partial

# Usage: python compute_dmps_pe.py [-p=num_proc] [-o=out_id] [-v=min_cov] [-c=chrm_list] <allc_path> <sample1> <sample2> [sampleN]*
NUMPROC=1
CHRMLIST=['Chr1','Chr2','Chr3','Chr4','Chr5']

def processInputs( allcPath, sampleNamesAr, outID, numProc, minCov, chrmList ):
	print( 'AllC path:', allcPath )
	print( 'Samples:', ', '.join(sampleNamesAr) )
	info = '#from_script: compute_dmps_pe.py'
	if minCov != None:
		info += '; min_cov: {:d}'.format( minCov )
	info += '; chrm_list: {:s}'.format(','.join(chrmList))
	info += '\n'
	filterDict = getLabelFilter( sampleNamesAr )
	
	# loop through chromosomes
	pool = multiprocessing.Pool( processes=numProc )
	
	results = [ pool.apply_async(processChrm, args=(allcPath, sampleNamesAr, filterDict, chrm, minCov) ) for chrm in chrmList ]
	dfAr = [p.get() for p in results]
	dfFish = pd.concat( dfAr, ignore_index=True )
	
	# correct p-values
	print( 'Correcting p-values ({:d} tests)'.format( dfFish.shape[0] ) )
	dfFish['pvalue.adjust'] = pValueAdjust( dfFish['pvalue'] )
	
	# reorder columns
	dfOut = dfFish[ ['sample_y','chrm', 'pos', 'mCounts_x', 'tCounts_x', 'mCounts_y', 'tCounts_y', 'pvalue', 'pvalue.adjust']]
	
	if outID == None:
		outID = 'out'
	outFileStr = outID + '_fisher_dmps.tsv'
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info)
	dfOut.to_csv( outFileStr, sep='\t', mode='a' )
	
	print( 'Done' )

def processChrm( allcPath, sampleNamesAr, filterDict, chrm, minCov ):
	print( 'Beginning', chrm )
	dfCat = pd.DataFrame()
	if os.path.exists('tmp_'+chrm+'_filter.csv')==False:
		# get allc files and concat
		for sampleName in sampleNamesAr:
			df = readAllc( allcPath, sampleName, chrm, minCov )
			dfCat = pd.concat( [dfCat, df], ignore_index=True )
		# filter for methylated positions
		print( 'Filtering',chrm,'for positions methylated in 1+ samples' )
		dfCat.index = dfCat['sample']
		dfGroup = dfCat.groupby( 'pos' )
		dfFilter = dfGroup.filter( lambda x: x['isMeth'].sum() > 0 )
		dfFilter.to_csv('tmp_'+chrm+'_filter.csv')
		del dfCat, dfGroup
	else:
		dfFilter = pd.read_csv('tmp_'+chrm+'_filter.csv')
		dfFilter.set_index('sample', inplace=True )
		dfFilter.rename( columns={'sample.1':'sample'}, inplace=True )
	
	# group by sample
	dfSample = dfFilter.groupby('sample')
	
	# call outer loop function
	print( 'Analyzing {:s}'.format( chrm ) )
	mapfunc = partial( runFisher, df=dfFilter, flt=filterDict )
	dfFish = dfSample.apply( runFisher, df = dfFilter, flt = filterDict )
	#res = pool.map( mapfunc, [g for g in dfSample] )
	#dfFish = pd.concat( res )
	dfFish['chrm'] = chrm
	dfFish.set_index( 'sample_x', inplace=True )
	del dfSample
	return dfFish

def readAllc( allcPath, sampleName, chrm, minCov ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}{:s}_{:s}.tsv'.format( allcPath, sampleName , ('' if minCov == None else '_cov'+str(minCov)), chrm ) )
	print( '~{:s}'.format(os.path.basename(inFileStr)))
	df = pd.read_table( inFileStr,header=None,names=['chrm', 'pos', 'strand', 'mcClass','mCounts', 'tCounts', 'isMeth'],comment='#',usecols=['pos','mCounts','tCounts','isMeth'] )
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

def runFisher( grDf, df, flt ):
	curGroup = grDf.index.unique()[0]
	#sample = group[0]
	#grDf = group[1]
	fltList = flt[curGroup]
	dfMerge = pd.merge( grDf, df, on='pos')
	# get only comparisons that we want and drop unnecessary columns
	dfFilter = dfMerge.query( 'sample_y in @fltList' )
	dfFilterD = dfFilter.drop( ['sample_x','isMeth_x','isMeth_y'],axis=1 )
	del dfMerge, dfFilter
	# group and apply
	dfComp = dfFilterD.groupby( 'pos'  )
	dfOutput = dfComp.apply( fisher )
	dfOutput['sample_x'] = curGroup
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
