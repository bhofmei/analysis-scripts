import sys, math, glob, multiprocessing, subprocess, os
import pandas as pd
import numpy as np
import scipy.stats as stats
from functools import partial
import pandas_util

# Usage: python compute_dmps_hdf5.py [-p=num_proc] [-o=out_id] [-v=min_cov] [-c=chrm_list] <allc_path> <sample1> <sample2> [sampleN]*
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
	
	if outID == None:
		outID = 'out'
		
	# check for hdf5 file
	hdfFileStr = 'store_'+outID+'.h5'
	store = pd.HDFStore( hdfFileStr )
	# does is fishers table there?
	if pandas_util.hdfTable( store, search='fisher', isObj=True) == False:
		# read all of the allc info
		readAllc( store, allcPath, sampleNamesAr, minCov )
		
		# filter allc
		processFilterMeth( store, chrmList )
		
		# outer fisher loop
		processFisher( store, chrmList, filterDict, sampleNamesAr )
	
	dfFish = store.select( 'fisher' )
	dfFish.set_index( 'sample_x', inplace=True )
	store.close()
	# correct p-values
	print( 'Correcting p-values ({:d} tests)'.format( dfFish.shape[0] ) )
	dfFish['pvalue.adjust'] = pandas_util.pvalueAdjust( dfFish['pvalue'] )
	
	# reorder columns
	dfOut = dfFish[ ['sample_y','chrm_x', 'pos', 'mCounts_x', 'tCounts_x', 'mCounts_y', 'tCounts_y', 'pvalue', 'pvalue.adjust']]
	
	store.close()
	outFileStr = outID + '_fisher_dmps.tsv'
	print( 'Writing output to', outFileStr )
	with open( outFileStr, 'w' ) as f:
		f.write(info)
	dfOut.to_csv( outFileStr, sep='\t', mode='a' )
	
	print( 'Done' )

def readAllc( storeObj, allcPath, sampleNameAr, minCov ):
	# check if allc table exists
	tableName = 'allcmulti'
	dfName = 'dfmulti'
	if pandas_util.hdfTable( storeObj, isObj=True, search=tableName):
		print( 'Reading allc from HDF file' )
		return True
	
	# first sample name need to initialize table
	df = readAllcFile( allcPath, sampleNameAr[0], minCov )
	storeObj.append_to_multiple({tableName:['chrm','pos','sample'], dfName: None}, df, selector=tableName)
	# loop through remaining
	for i in range(1, len(sampleNameAr) ):
		sampleName = sampleNameAr[i]
		df = readAllcFile( allcPath, sampleName, minCov )
		storeObj.append_to_multiple({tableName:['chrm','pos','sample'], dfName: None}, df, selector=tableName)
	# end for i
	
def readAllcFile( allcPath, sampleName, minCov ):
	inFileStr = os.path.normpath( '{:s}/allc_{:s}{:s}.tsv'.format( allcPath, sampleName , ('' if minCov == None else '_cov'+str(minCov)) ) )
	print( ' {:s}'.format(os.path.basename(inFileStr)))
	df = pd.read_table( inFileStr,header=None,names=['chrm', 'pos', 'strand', 'mcClass','mCounts', 'tCounts', 'isMeth'],comment='#',usecols=['chrm','pos','mCounts','tCounts','isMeth'] )
	df['sample'] = sampleName
	return df
	
def processFilterMeth( storeObj, chrmList ):
	# loop through chromosomes
	tableName = 'allcmultifilter'
	dfName = 'dfmultifilter'
	tableName2 = 'allcmulti'
	dfName2 = 'dfmulti'
	if pandas_util.hdfTable( storeObj, isObj=True, search=tableName):
		print( 'Reading filtered allc from HDF file' )
		return True
	# else do the filtering
	for cc in chrmList:
		print( 'Filtering', cc )
		dfP = storeObj.select_as_multiple( [tableName2,dfName2], where=["chrm=cc"], selector=tableName2 )
		dfG = dfP.groupby( 'pos' )
		dfF = dfG.filter( lambda x: x['isMeth'].sum() > 0 )
		storeObj.append_to_multiple({tableName:['chrm','pos','sample'], dfName: None}, dfF, selector=tableName)
		# list of positions
		'''posList = storeObj.select( tableName2, where="chrm=cc", columns=['pos'] )['pos'].unique()
		n = len(posList)
		nn = n // 10
		for i in range(n):
			if posList[i] % nn == 0:
				print( '{:s}%..'.format( posList[i] // nn ), end='' )
			spo=str(posList[i])
			# get data, pass to function, add if necessary
			df = storeObj.select_as_multiple( [tableName2,dfName2], where=["chrm=cc","pos=spo"], selector=tableName2 )
			if df['isMeth'].sum() > 0:
				storeObj.append_to_multiple({tableName:['chrm','pos','sample'], dfName: None}, df, selector=tableName)
		# end for position
		'''
		print()
	# end for chrm

def processFisher( storeObj, chrmList, filterDict, sampleNamesAr ):
	tableName = 'fisher'
	if pandas_util.hdfTable( storeObj, isObj=True, search=tableName):
		print( 'Reading fisher test results from HDF file' )
		return True
	
	# loop through chromosome
	for cc in chrmList:
		print( 'Analyzing', cc )
		# loop through samples
		for sampleName in sampleNamesAr:
			print( '~', sampleName )
			processSample( storeObj, cc, sampleName, filterDict[sampleName] )
		# end for sample
	#end for cc
	# reindex
	storeObj.create_table_index( tableName, columns=['sample_x'] )

def processSample( storeObj, cc, sampleName, compNameAr ):
	tableName1 = 'allcmultifilter'
	dfName = 'dfmultifilter'
	tableName2 = 'fisher'
	# this sample's data
	df = storeObj.select_as_multiple([tableName1, dfName], where=["chrm=cc",'sample=sampleName'], selector=tableName1)
	for flt in compNameAr:
		print( '~~', flt )
		dff = storeObj.select_as_multiple([tableName1, dfName], where=["chrm=cc",'sample=flt'], selector=tableName1)
		dfMerge = pd.merge( df, dff, on='pos')
		dfMerge.drop( ['isMeth_x', 'isMeth_y','chrm_y'], axis=1, inplace=True )
		dfComp = dfMerge.groupby( 'pos' )
		gfOut = dfComp.apply( fisher )
		storeObj.append( tableName2, gfOut, data_columns=['chrm','pos'], index=False )
	# end for flt
		
def getLabelFilter( labels ):
	filterDict = {}
	for i in range(len(labels)):
		filterDict[labels[i]] = []
		for j in range(i+1,len(labels)):
			filterDict[labels[i]] += [ labels[j] ]
		# end for j
	# end for i
	return filterDict

def fisher( data ):
	''' excutes fisher's exact test on incoming data
		returns original input data with additional column containing test
		p-value
	'''
	mCounts = np.asarray(data.loc[:,['mCounts_x','mCounts_y']])[0]
	tCounts = np.asarray(data.loc[:,['tCounts_x','tCounts_y']])[0]
	uCounts = tCounts - mCounts
	subData = np.array([mCounts,uCounts])
	odr, p = stats.fisher_exact( subData )
	data['pvalue']=p
	return data

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
