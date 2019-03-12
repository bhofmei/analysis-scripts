import pandas as pd
import numpy as np
import math
from decodingpath import *
from transitions import Transitions

EPS=0.0001

def runFB( group, trans ):
		fb = DecodeForwardBackward( group, trans )
		fbr = fb.run()
		#print( fbr.head() )
		#return pd.Series(bin=fbr['bin'],prediction=fbr['fb.prediction'])
		fbr['prediction'] = fbr['fb.prediction']
		return fbr.drop( ['fb.prediction', 'fb.score.mother', 'fb.score.MPV', 'fb.score.father'], axis=1)

class AdaptiveTransitionsv2:
	
	def __init__( self, df, trans, ignoreAr, iter ):
		'''smpName = emis.index.levels[0][emis.index.labels[0][0]]
		print( '  BW transitions for', smpName, end='...' )
		if smpName in ignoreAr:
			iter = 0'''
		# emis needs index reset
		self.labels = np.array( ['mother', 'MPV', 'father'] )
		self.data = df
		self.transitions = trans
		self.emissions = self._getEmissions( df )
		self.maxIter = iter
		self.states, self.size = self.emissions.shape
	
	def _getEmissions( self, data ):
		idVars = ['sample','bin','prediction']
		if 'num.feat' in list(data):
			idVars += ['num.feat']
		dfm = pd.melt( data, id_vars=idVars)
		dfp = dfm.pivot( index='variable', columns='bin', values='value' )
		dfe = dfp.reindex( self.labels )
		return dfe.values
		
	def run( self ):
		### going to try running forward-backward, updating transitions from
		# that until convergence
		#nSamples = self.observations - len(self.ignore)
		curIter = 0
		while curIter < self.maxIter:
			curIter += 1
			oldTransitions = np.copy( self.transitions )
			# get predictions
			#newPredictions = self.groups.apply( runFB, oldTransitions )
			newPredictions = runFB(self.data, oldTransitions )
			#print( newPredictions.head() )
			# get transitions
			t = Transitions( newPredictions, [] )
			self.transitions = t.getTransitions()
			#print( self.transitions )
			# check if transitions converged
			if np.linalg.norm( oldTransitions - self.transitions ) < EPS:
				break
		print( 'did not converge' if curIter == self.maxIter else 'coverged in {:d} iterations'.format( curIter ) )
		return pd.DataFrame(self.transitions,index=self.labels,columns=self.labels)


