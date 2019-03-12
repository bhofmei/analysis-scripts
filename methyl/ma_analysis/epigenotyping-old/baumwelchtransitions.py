import pandas as pd
import numpy as np
import math

EPS=0.00001

class TransitionBaumWelch:
	
	def __init__( self, emis, transMat, iter=10, ignoreAr=None ):
		# this is for a single sample with an initial transition matrix
		# need to get sample name from it first
		smpName = emis.index.levels[0][emis.index.labels[0][0]]
		print( '  BW transitions for', smpName, end='...' )
		if smpName in ignoreAr:
			iter = 0
		# emis needs index reset
		self.labels = np.array( ['mother', 'MPV', 'father'] )
		self.emissions = emis.reset_index( level=0,drop=True ).reindex( self.labels ).values
		self.emissions = np.exp( self. emissions )
		self.transitions = transMat
		self.maxIter = iter
		self.states, self.size = self.emissions.shape
	
	def _initialize( self ):
		self.theta = np.zeros( (self.states, self.states, self.size) )
		self.forward = np.zeros( (self.states, self.size+1) ) # m x n+1
		self.forward[:,0] = 1.0/self.states
		self.backward = np.zeros( (self.states, self.size+1) ) # m x n+1
		self.backward[:,-1] = 1.0
		# initialize posterior prob dist
		self.posterior = np.zeros( (self.size, self.states) ) # n x m
	
	def run( self ):
		self._initialize()
		curIter = 0
		while curIter < self.maxIter:
			curIter += 1
			print( curIter, self.transitions )
			oldTransitions = self.transitions
			newTransitions = self._runIter()
			self.transitions = np.copy( newTransitions )
			if np.linalg.norm(oldTransitions - self.transitions ) < EPS:
				break
		# end while
		print( 'did not converge' if curIter == self.maxIter else 'coverged in {:d} iterations'.format( curIter ) )
		return pd.DataFrame( self.transitions, index=self.labels, columns=self.labels )
		
	def _runIter( self ):
		# run forward-backward
		self._forwardBackward()
		
		# fill theta
		for i in range(self.states):
			for j in range(self.states):
				for k in range(self.size):
					self.theta[i,j,k] = self.forward[i,k] * self.backward[j,k+1] * self.transitions[i,j] * self.emissions[j,k]
				# end for k
			# end for j
		# end for i
		
		# new transition matrix
		newTrans = np.zeros( (self.states, self.states) )
		for i in range(self.states):
			for j in range(self.states):
				newTrans[i,j] = np.sum( self.theta[i,j,:]) / np.sum(self.posterior[i:])
			# end for j
		# end for i
		newTrans = newTrans / np.sum(newTrans, 1 )
		return newTrans
	
	def _forwardBackward( self ):
		# intialize first rows
		self.forward[:,0] = 1.0/self.states
		self.backward[:,-1] = 1.0
		# fill forward -> loop across bins
		for i in range(self.size):
			# get current column values
			fCol = np.matrix( self.forward[:,i] )
			# fill in next column
			self.forward[:,i+1] = fCol * np.matrix( self.transitions ) * np.matrix( np.diag( self.emissions[:,i] ) )
			# normalize
			self.forward[:,i+1] = self.forward[:,i+1] / np.sum( self.forward[:,i+1] )
		# end for i
		
		# fill backwards -> loop across bins
		for i in range( self.size, 0, -1 ):
			# get current column values
			bRow = np.matrix( self.backward[:,i]).transpose()
			# get values for next column
			tmpCol = ( np.matrix(self.transitions) * np.matrix(np.diag(self.emissions[:,i-1])) * bRow).transpose()
			# normalize
			self.backward[:,i-1] = tmpCol / np.sum( tmpCol )
		# end for i
		
		# combine
		tmpPosterior = np.zeros((self.states, self.size))
		tmpPosterior = np.array( self.forward ) * np.array( self.backward )
		#tmpPosterior = np.array( self.forward[:,1:] ) * np.array( self.backward[:,:-1] )
		# normalize
		tmpPosterior = tmpPosterior / np.sum( tmpPosterior, 0)
		self.posterior = np.transpose(tmpPosterior)
