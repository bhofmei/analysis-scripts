import numpy as np
import pandas as pd
import math

# module for computing transition probabilities
# simple transitions independent of time in each state

class SimpleTransitions:
	
	def __init__( self, df, ignore=[] ):
		self.data = df
		self.transCountDict = {}
		self.labels = ['mother', 'MPV', 'father', 'total']
		self.ignoreLabels = ignore
		self.transProbMat = []
	
	def run( self ):
		self._initialize()
		self._fillDict()
		self._computeProbs()
		return self.transProbMat
		
		
	def _initialize( self ):
		self.transCountDict = {}
		for x in self.labels:
			self.transCountDict[x] = {}
			for y in self.labels:
				self.transCountDict[x][y] = 1
			# end for y
		# end for x
		self.transProbMat = []
		
	def _fillDict( self ):
		dfp = self.data.pivot( index='sample', columns = 'bin', values='prediction' )
		for index,row in dfp.iterrows():
			if index not in self.ignoreLabels:
				for i in range(row.size - 1):
					try:
						tmpDict = self.transCountDict[row.iloc[i]]
						tmpDict[row.iloc[i]] += 1
						tmpDict['total'] += 1
					except KeyError:
						print( 'key', row.iloc[i], 'or', row.iloc[i+1], 'not found' )
				# end for i
		# end for index, row
	
	def _computeProbs( self ):
		for x in self.labels:
			ar = [ math.log(float(self.transCountDict[x][y]) / self.transCountDict[x]['total']) for y in self.labels ]
			self.transProbMat.append( ar )
	
	def __str__( self ):
		if self.transProbMat == []:
			return 'transitions not computed'
		outStr = ''
		lAr = ['m','h','f']
		headerAr = [ ' {:>9s}'.format( x ) for x in lAr ]
		outStr += ' '*2 + ''.join(headerAr) + '\n' + ' '*3 + '-'*(3*10) +'\n'
		for i in range(len(lAr)):
			tmpStr = ' {:s}|'.format( lAr[i] )
			transAr = [ ' {:>9.6f}'.format(self.transProbMat[i][j]) for j in range(3) ]
			outStr += tmpStr + ''.join(transAr)+ '|\n'
		outStr += ' '*3 + '-'*(3*10)
		return outStr
