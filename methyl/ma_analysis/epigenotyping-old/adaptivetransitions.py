import pandas as pd
import numpy as np
import math
from decodingpath import *
from transitions import Transitions

EPS=0.0001

def runFB( group, trans ):
	print('-\n', group.shape)
	fb = DecodeForwardBackward( group, trans )
	fbr = fb.run()
	#print( fbr.head() )
	#return pd.Series(bin=fbr['bin'],prediction=fbr['fb.prediction'])
	fbr['prediction'] = fbr['fb.prediction']
	return fbr.drop( ['fb.prediction', 'fb.score.mother', 'fb.score.MPV', 'fb.score.father'], axis=1)

class AdaptiveTransitions:
	
	def __init__( self, df, trans, ignoreAr, iter ):
		# n bins, m states
		self.labels = np.array( ['mother', 'MPV', 'father'] )
		self.data = pd.DataFrame( df.values.copy(), df.index.copy(), df.columns.copy() )
		self.transitions = trans # fraction probs not log, m x m
		self.groups = self.data.groupby( 'sample' )
		self.observations = len( self.groups.groups )
		#print(self.observations)
		self.states = len( self.labels )
		self.size = self.data.shape[0] / self.observations
		self.maxIter = iter
		self.ignore = ignoreAr

	def run( self ):
		### going to try running forward-backward, updating transitions from
		# that until convergence
		#nSamples = self.observations - len(self.ignore)
		curIter = 0
		while curIter < self.maxIter:
			curIter += 1
			oldTransitions = np.copy( self.transitions )
			# get predictions
			newPredictions = self.groups.apply( runFB, oldTransitions )
			#print( newPredictions.head() )
			# get transitions
			t = Transitions( newPredictions, self.ignore )
			self.transitions = t.getTransitions()
			#print( self.transitions )
			# check if transitions converged
			if np.linalg.norm( oldTransitions - self.transitions ) < EPS:
				break
		return curIter, self.transitions
	
		
