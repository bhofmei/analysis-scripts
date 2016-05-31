import pandas as pd
import numpy as np

class ProbPath:
	
	def __init__(self, df, transitions ):
		self.transMat = transitions
		self.data = df
		self.size = self.data.shape[0]
		self.probMat = None
		self.traceMat = None
		self._initializeMat()
	
	def _initializeMat( self ):
		nbin = self.size
		self.probMat = [ [0]*3 for x in range(nbin) ]
		self.traceMat = [ [None]*3 for x in range(nbin) ]
	
	def run( self ):
		self._fillMat()
		#self._printMat()
		self._followTraceback()
		#self._printTraceMat()
		return self.data
	
	def _fillMat( self ):
		labels = ['mother','MDV','father']
		# loop through rows/bin
		for i in range(self.size):
			# loop through columns/states
			for j in range(3):
				prob = self.data[labels[j]].iloc[i]
				maxS, maxP = self._computeScore( i, j, prob )
				#print( prob, maxS )
				self.probMat[i][j] = maxS
				self.traceMat[i][j] = maxP
			# end for j
		# end for i
	
	def _computeScore( self, i, j, prob ):
		scores = [prob]*3
		for k in range(3):
			#scores[k] *= ( 1 if i == 0 else self.probMat[i-1][k]) * self.transMat[k][j]
			# change this to use log scores
			scores[k] += ( 0 if i == 0 else self.probMat[i-1][k]) + self.transMat[k][j]
		#print( scores )
		maxS = max( scores )
		maxI = scores.index( maxS )
		maxPos = ( None if i == 0 else maxI )
		return maxS, maxPos
	
	def _printMat( self ):
		outStr = ''
		lAr = ['m','h','f']
		headerAr = [ ' {:>8s}'.format( x ) for x in lAr ]
		outStr += ' '*5 + ''.join(headerAr) + '\n' + ' '*5 + '-'*(3*9) +'\n'
		for i in range(20):
			tmpStr = '{:<4d}|'.format( i )
			transAr = [ ' {:>8.1f}'.format(self.probMat[i][j]) for j in range(3) ]
			outStr += tmpStr + ''.join(transAr)+ '|\n'
		outStr += ' '*5 + '-'*(3*9)
		print( outStr )
	
	def _followTraceback( self ):
		# add columns to dataframe for output
		self.data['score'] = 0
		self.data['iprediction'] = 'NA'
		vals = self.probMat[self.size-1]
		maxS = max( vals )
		nextJ = vals.index( maxS )
		for i in range(self.size-1, -1, -1):
			nextJ = self._tracebackHelper( i, nextJ )
			if nextJ == None:
				break
		# end for i
	
	def _tracebackHelper( self, i, j ):
		sj = np.nonzero(self.data.columns.values == 'score')[0][0]
		si = np.nonzero(self.data.columns.values == 'iprediction')[0][0]
		labels = ['mother','MDV','father']
		score = self.probMat[i][j]
		label = labels[j]
		self.data.iloc[i,sj] = score
		self.data.iloc[i,si] = label
		return self.traceMat[i][j]
	
	def _printTraceMat( self ):
		outStr = ''
		lAr = ['m','h','f']
		lAr2 = ['mother', 'MPV', 'father' ]
		headerAr = [ ' {:>6s}'.format( x ) for x in lAr ]
		outStr += ' '*5 + ''.join(headerAr) + '\n' + ' '*5 + '-'*(3*7) +'\n'
		for i in range(20):
			tmpStr = '{:<4d}|'.format( i )
			transAr = [ ' {:>6s}'.format(lAr2[self.traceMat[i][j][1]]) for j in range(3) ]
			outStr += tmpStr + ''.join(transAr)+ '|\n'
		outStr += ' '*5 + '-'*(3*7)
		print( outStr )
		
