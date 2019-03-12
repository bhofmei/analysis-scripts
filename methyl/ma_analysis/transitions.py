### Transition matrix

import numpy as np
import pandas as pd
import math

class Transitions:
	'''
		compute transition probabilities
	'''
	def __init__( self, df, ignore=[] ):
		self.data = df
		self.labels = ['mother', 'MPV', 'father', 'total' ]
		self.ignore = ignore
	
	def getTransitions( self ):
		''' function to be access by outside classes '''
		self._initialize()
		self._fillDict()
		self._computeProbs()
		return np.array( self.transMatrix )[0:3,0:3]
	
	def _initialize( self ):
		''' initialize the counting dictionary '''
		self.countDict = {}
		for x in self.labels:
			self.countDict[x] = {}
			for y in self.labels[:-1]:
				self.countDict[x][y] = 1
			self.countDict[x]['total'] = len(self.labels) - 1
			# end for y
		# end for x
		self.transMatrix = []
	
	def _fillDict( self ):
		''' fill count dictionary '''
		dfp = self.data.pivot( index='sample', columns = 'bin', values = 'prediction' )
		for index,row in dfp.iterrows():
			#print( index, row[:5] )
			if index not in self.ignore:
				for i in range(row.size-1):
					try:
						tmpDict = self.countDict[row.iloc[i]]
						tmpDict[row.iloc[i+1]] += 1
						tmpDict['total'] += 1
					except KeyError:
						print( 'key', row.iloc[i], 'or', row.iloc[i+1], 'not found' )
				# end for i
		# end for index,row
	
	def _computeProbs( self ):
		for x in self.labels:
			ar = [ float(self.countDict[x][y]) / self.countDict[x]['total'] for y in self.labels ]
			self.transMatrix.append( ar )
