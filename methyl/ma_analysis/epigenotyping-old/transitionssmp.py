import pandas as pd
import numpy as np
import math
from decodingpath import *

EPS=0.000001

class ImproveTransitions:
	
	def __init__( self, data, transMat, iter ):
		self.labels = ['mother', 'MPV', 'father' ]
		self.data = pd.DataFrame( data.values.copy(), data.index.copy(), data.columns.copy() )
		self.transitions = transMat
		self.maxIter = iter
		self.size = self.data.shape[0]
		self.states = len(self.labels)

	def getTransitions( self ):
		curIter = 0
		while curIter < self.maxIter:
			curIter += 1
			oldTransitions = np.copy( self.transitions )
			#print( '-', self.transitions.dtype)
			fb = DecodeForwardBackward( self.data, oldTransitions )
			newPredictions = fb.run()
			#newPredictions['prediction'] = newPredictions['fb.prediction']
			pred = newPredictions['fb.prediction'].values
			#print(pred)
			# update transition matrix
			self.transitions = self._newTransitionMatrix( pred )
			if np.linalg.norm(oldTransitions - self.transitions ) < EPS:
				break
		# end while
		#print( self. transitions )
		return self.transitions, curIter
			
	def _newTransitionMatrix( self, pred ):
		# initialize array
		#print(self.size)
		countAr = np.ones( (self.states,self.states) )
		for i in range(self.size-1):
			l= pred[i:i+2]
			j,k = [ self.labels.index( x ) for x in (l) ]
			#print('*',j,k)
			countAr[j,k] += 1
		# iterate other direction
		for i in range(self.size-1, 0, -1):
			l = [ pred[i] ] + [ pred[i-1] ]
			j,k = [ self.labels.index( x ) for x in (l) ]
			#print('**',j,k)
			countAr[j,k] += 1
		# fix proportions
		#print(countAr)
		return countAr / np.sum(countAr, 1, keepdims=True )
