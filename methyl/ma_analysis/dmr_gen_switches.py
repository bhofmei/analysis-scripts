import sys, math, glob, multiprocessing, subprocess, os, bisect, random
import pandas as pd
import numpy as np
import scipy.stats as stats
from functools import partial

NUMC=10
METH=0.3
FDR=0.05
LENT=40

# Usage: python dmr_gen_switches.py [-wm] [-n=num_c_thresh] [-m=meth_thresh] [-f=fdr] [-d=length_thresh] [-o=outID] <in_file>

def processInputs(inFileStr, outID, numC, mThresh, fdrThresh, lenThresh, isRawWM ):
	print( 'Input file', os.path.basename( inFileStr ) )
	print( 'Number of C\'s threshold:', numC )
	
	if isRawWM:
		print( 'Weighted methylation threshold:', mThresh )
	else:
		print( 'Methylation change threshold:', mThresh )
	print( 'Length threshold:', lenThresh )
	print( 'FDR:', fdrThresh )
	
	# read input file into a table
	print( ' Reading', os.path.basename(inFileStr) )
	df = pd.read_csv(inFileStr, comment='#',sep='\t')
	df.set_index('DMR',inplace=True)
	print( ' Found {:d} DMRs'.format(len(df.index.unique()) ) )
	print( ' Filtering' )
	dfD = df.groupby( df.index )
	dfF = dfD.filter( filterDMR, nc=numC, wm=mThresh, ln=lenThresh, ir=isRawWM )
	dfFC = dfF.copy()
	del df, dfD, dfF
	
	nDMRs = len(dfFC.index.unique())
	print( ' Analyzing remaining {:d} DMRs'.format( nDMRs ) )
	dfFC['pvalue'] = dfFC.apply( computePropTest, axis=1 )
	
	print( 'Applying correction for {:d} tests'.format( dfFC.shape[0] ) )
	dfFC['pvalue.adjust'] = pvalueAdjust( dfFC['pvalue'] )
	dfFC['pvalue.log'] = -1 * np.log10( dfFC['pvalue.adjust'] )
	
	print( 'Analyzing for switches' )
	dfFC['wm.change'] = dfFC.apply( lambda x: ( float(x['wm1'] - x['wm2']) / float( (1.0 if x['wm1']==0 else x['wm1']) ) ), axis=1 )
	if isRawWM:
		dfFC['mThresh'] = mThreshold(dfFC['wm.diff'],wm=mThresh)
	else:
		dfFC['mThresh'] = mThreshold( dfFC['wm.change'], wm=mThresh )
	logThresh = -1 * math.log10( fdrThresh )
	dfFC['pThresh'] = pThreshold( dfFC['pvalue.log'], pt=logThresh )
	dfFC['isSwitch'] = dfFC.apply( lambda x: (x['mThresh'] and x['pThresh']), axis=1)
	
	print( 'Counting switches' )
	dfFG = dfFC.groupby( [dfFC.index, 'region'] )
	dfSw = dfFG.agg( {'isSwitch': np.sum} )
	dfSw.rename(columns={'isSwitch':'switch.counts'},inplace=True)
	
	dfSF = pd.DataFrame( {'freq': np.bincount( dfSw['switch.counts'] ) } )
	dfSF.index.name = 'switch.counts'
	
	# write outputs
	dfSw.reset_index( inplace=True )
	if outID == None:
		t = os.path.basename( inFileStr )
		outID = t.replace('.tsv','')
	info = '#from_script: dmr_gen_switches.py; in_file: {:s}; num_c_thresh: {:d}; wei_meth_thresh: {:g}; raw_wei_meth: {:s}; fdr: {:g}; len_thresh: {:d}; dmr_count: {:d}\n'.format( os.path.basename( inFileStr ), numC, mThresh, str(isRawWM), fdrThresh, lenThresh, nDMRs )
	outFileStrAr = [ outID + '_' + x + '.tsv' for x in ['full', 'switches', 'switch_counts' ] ]
	print( 'Writing outputs to', ', '.join( outFileStrAr ) )
	for outFileStr in outFileStrAr:
		with open( outFileStr, 'w' ) as f:
			f.write(info)
	dfFC.to_csv( outFileStrAr[0], sep='\t', mode='a' )
	dfSw.to_csv( outFileStrAr[1], sep='\t', mode='a' )
	dfSF.to_csv( outFileStrAr[2], sep='\t', mode='a' )
	print( 'Done' )
	
def filterDMR( data, nc, wm, ln, ir ):
	''' nc = number of cytosines
		wm = weighted methylation
		ir = is raw methylation
	'''
	if data['num.cs'].min() < nc:
		return False
	if data['length'].min() < ln:
		return False
	#r = np.append( data['mC.r1'], data['mC.r2'] )
	#if r.min() < 1:
	#	return False
	w = np.append( data['wm1'], data['wm2'] )
	a = w.max() - w.min()
	if not ir:
		a = a / w.min()
	if a < wm:
		return False
	return True

def computePropTest( data ):
	mC = np.asarray( data[['mC.r1', 'mC.r2']] )
	tC = np.asarray( data[['t.r1','t.r2']] )
	uC = tC - mC
	#chi2, p, dof, ex = stats.chi2_contingency( np.array([mC,uC]), correction=True)
	odr, p = stats.fisher_exact( np.array([mC,uC]) )
	return p

def pvalueAdjust( pvalues ):
	''' expects a list/series/array of p-values
		returns list of adjust pvalues
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

def mThreshold( data, wm=0.3, ir=False ):
	n = len(data)
	out = np.zeros(n, dtype=np.int)
	for i, m in enumerate(data):
		out[i] = ( abs(m) >= wm )
	return out

def pThreshold( data, pt ):
	n = len(data)
	out = np.zeros(n, dtype=np.int)
	for i, m in enumerate(data):
		out[i] = m >= pt
	return out

def parseInputs( argv ):
	numC = NUMC
	mThresh = METH
	fdrThresh = FDR
	lenThresh = LENT
	isRawWM= False
	outID = None
	startInd = 0
	
	for i in range(min(4,len(argv))):
		if argv[i].startswith( '-o=' ):
			outID = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-m=' ):
			try:
				mThresh = float( argv[i][3:] )
				startInd += 1
				if mThresh > 1:
					mThresh = mThresh / 100.0
			except ValueError:
				print( 'WARNING: methylation threshold must be numeric...using', METH )
				mThresh = METH
		elif argv[i] == '-wm':
			isRawWM = True
			startInd += 1
		elif argv[i].startswith( '-f=' ):
			try:
				fdrThresh = float( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: FDR must be numeric...using', FDR )
				fdrThresh = FDR
		elif argv[i].startswith( '-d=' ):
			try:
				lenThresh = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'WARNING: length threshold must be integer...using', LENT )
				lenThresh = LENT
	# end for
	inFileStr = argv[startInd]
	processInputs(inFileStr, outID, numC, mThresh, fdrThresh, lenThresh, isRawWM )

def printHelp():
	print( 'Usage: python dmr_gen_switches.py [-wm] [-n=num_c_thresh] [-m=meth_thresh] [-d=length_thresh] [-f=fdr] [-o=outID] <in_file>' )
	print( 'Required:' )
	print( 'in_file\ttab-delimited file of DMRs and read counts' )
	print( 'Optional:' )
	print( '-wm\tmethylation threshold is for raw methyl difference\n\t\tnot percent difference' )
	print( '-n=num_c_thresh\tmin number of cytosines in region to be considered\n\t\tfor analysis [default 10]' )
	print( '-d=lenth_thresh\tmin length of dmr in bp [default 40]' )
	print( '-m=meth_thresh\tmin methylation change btwn generations to be\n\t\tconsidered a switch [default 0.3]' )
	print( '-f=fdr\tFDR value for significant switches [default 0.05]' )
	print( '-o=out_id\tidentifier for output files [default uses input file name]' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
