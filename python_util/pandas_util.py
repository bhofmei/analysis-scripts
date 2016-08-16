# pandas util functions
import numpy as np

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

def hdfTable( inFile, search=None, isObj=False ):
	'''
		when search not specified, returns list of table names
		when search is specified, returns True is search is in file
	'''
	if isObj:
		k = inFile.keys()
	else:
		k = pd.HDFStore( inFile ).keys()
	ks = [ ( x[1:] if x.startswith( '/' ) else x ) for x in k ]
	if search == None:
		return ks
	return ( search in ks )

def hdfTableList( inFile, searchList, isObj ):
	if isObj:
		k = inFile.keys()
	else:
		k = pd.HDFStore( inFile ).keys()
	ks = [ ( x[1:] if x.startswith( '/' ) else x ) for x in k ]
	for sl in searchList:
		if sl in ks:
			return True
	# end for sl
	return False
