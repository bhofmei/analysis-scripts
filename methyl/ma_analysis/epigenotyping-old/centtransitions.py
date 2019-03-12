import pandas as pd
import numpy as np
import math

class Transitions:
	
	def __init__( self, df, ignore=[], cent=[] ):
		self.data = df
		self.labels = ['mother', 'MPV', 'father']
		self.ignore = ignore
		self.centromere=cent
		
		# initialize count matrix
		self.transitions = np.ones( (len(self.labels), len(self.labels)) )
	
	def _fill( self ):
		dfp = self.data.pivot( index='sample', columns='bin', values='prediction' )
		colnames = list( dfp )
		for index, row in dfp.iterrows():
			if index not in self.ignore:
				for i in range(len(colnames)-1):
					# ignore those in centromere
					if colnames[i] in self.centromere or colnames[i+1] in self.centromere:
						continue
					# get state index
					l1 = self.labels.index( row[colnames[i]] )
					l2 = self.labels.index( row[ colnames[i+1] ] )
					self.transitions[l1,l2] += 1
				# end for i
		# end for index
		
		# normalize
		self.transitions = self.transitions / np.sum(self.transitions, 1, keepdims=True)
	
	def getTransitions( self ):	
		self._fill()
		return self.transitions
