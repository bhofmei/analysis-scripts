# Transitions class
import numpy as np
import pandas as pd

class Transitions:
	
	def __init__( self, df, ignore=[] ):
		self.data = df
		self.labels =  ['mother', 'MPV', 'father', 'total']
		self.ignoreLabels = ignore
	
	def getTransitions( self ):
		self._initialize()
		self._fillDict()
		self._computeProbs( )
		print( self.transCountDict )
		return np.array(self.transProbMat)[0:3,0:3]
	
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
						tmpDict[row.iloc[i+1]] += 1
						tmpDict['total'] += 1
					except KeyError:
						print( 'key', row.iloc[i], 'or', row.iloc[i+1], 'not found' )
				# end for i
		# end for index, row
	
	def _computeProbs( self ):
		for x in self.labels:
			ar = [ float(self.transCountDict[x][y]) / self.transCountDict[x]['total'] for y in self.labels ]
			self.transProbMat.append( ar )
	
	
