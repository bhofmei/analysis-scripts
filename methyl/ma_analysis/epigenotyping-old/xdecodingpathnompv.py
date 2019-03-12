### Decoding types ###
import pandas as pd
import numpy as np
import math

class DecodePath:
	''' Base class for decoding types
	'''
	def __init__( self, df, transMat ):
		# n bins, m states
		self.labels = np.array( ['mother', 'father'] )
		self.data = df	# n x m
		self.transitions = transMat # fraction probs not log, m x m
		self.emissions = self._getEmissions( df )	# m x n
		self.states, self.size = self.emissions.shape
	
	def _getEmissions( self, data ):
		dfm = pd.melt( data, id_vars=['sample','bin','prediction'])
		dfp = dfm.pivot( index='variable', columns='bin', values='value' )
		dfe = dfp.reindex( self.labels )
		return dfe.values

class DecodeViterbi( DecodePath ):
	''' Viterbi decoding
	'''
	def run( self ):
		''' main function accessible by outside classes '''
		self._initializeV()
		self._fillV()
		self._pathV()
		return self.data
	
	def _initializeV( self ):
		# take log of transitions
		self.log_transitions = np.log( self.transitions )	# m x m
		# initialize empty data structures for dynamic programming
		self.probabilities = np.zeros( (self.size, self.states) )	# n x m
		self.traceback = np.zeros( (self.size, self.states), dtype=np.int8 ) # n x m
	
	def _fillV( self ):
		# loop through rows/bins
		for i in range(self.size):
			# loop through states
			for j in range(self.states):
				em = self.emissions[j,i]	# note: m x n
				maxS, maxP = self._computeScore( i, j, em )
				self.probabilities[i,j] = maxS
				self.traceback[i,j] = maxP
			# end for j
		# end for i
	
	def _computeScore( self, i, j, prob ):
		scores = np.array( [prob]*self.states )	# 1 x m
		for k in range(self.states):
			if i != 0:
				scores[k] += self.probabilities[i-1,k]
			scores[k] += self.log_transitions[k,j]
		# end for k
		maxS = scores.max()
		maxP = (-1 if i == 0 else scores.argmax() )
		return maxS, maxP
	
	def _pathV( self ):
		# add columns to output
		self.data['vit.score.mother'] = self.probabilities[:,0]
		#self.data['vit.score.MPV'] = self.probabilities[:,1]
		self.data['vit.score.father'] = self.probabilities[:,1]
		self.data['vit.prediction'] = 'NA'
		vals = self.probabilities[self.size-1]
		# start traceback
		nextJ = vals.argmax()
		for i in range( self.size-1, -1, -1):
			nextJ = self._tracebackHelper( i, nextJ )
			if nextJ == -1:
				break # finished traceback
		# end for i
	
	def _tracebackHelper( self, i, j ):
		# get column numer where to record decoded prediction
		colI = np.nonzero(self.data.columns.values == 'vit.prediction')[0][0]
		# get current label to record
		label = self.labels[j]
		self.data.iloc[i, colI] = label
		# return next cell to travel to
		return self.traceback[i,j]
	

class DecodeForwardBackward( DecodePath ):
	''' Forward-backward decoding
	'''
	def run( self ):
		''' main function accessible by outside classes '''
		self._initializeF()
		self._fillF()
		self._pathF()
		return self.data
	
	def _initializeF( self ):
		# transform emissions from log to fractions
		self.prob_emissions = np.exp( self.emissions )
		# initialize forward and backward dynamic programming structures
		self.forward = np.zeros( (self.states, self.size+1) ) # m x n+1
		self.backward = np.zeros( (self.states, self.size+1) ) # m x n+1
		self.backward[:,-1] = 1.0
		# initialize posterior prob dist
		self.posterior = np.zeros( (self.size, self.states) ) # n x m
	
	def _fillF( self ):
		# fill forward -> loop across bins
		for i in range(self.size):
			# get current column values
			fCol = np.matrix( self.forward[:,i] )
			# fill in next column
			tmpCol = fCol * np.matrix( self.transitions) * np.matrix( np.diag( self.prob_emissions[:,i] ) )
			# normalize
			self.forward[:,i+1] = tmpCol / np.sum( tmpCol )
		# end for i
		
		# fill backwards -> loop across bins
		for i in range( self.size, 0, -1 ):
			# get current column values
			bRow = np.matrix( self.backward[:,i]).transpose()
			# get values for next column
			tmpCol = ( np.matrix(self.transitions) * np.matrix(np.diag(self.prob_emissions[:,i-1])) * bRow).transpose()
			# normalize
			self.backward[:,i-1] = tmpCol / np.sum( tmpCol )
		# end for i
		
		# combine
		tmpPosterior = np.zeros((self.states, self.size))
		tmpPosterior = np.array( self.forward[:,1:] ) * np.array( self.backward[:,:-1] )
		# normalize
		tmpPosterior = tmpPosterior / np.sum( tmpPosterior, 0)
		self.posterior = np.transpose(tmpPosterior)
	
	def _pathF( self ):
		# add columns to output
		self.data['fb.score.mother'] = self.posterior[:,0]
		#self.data['fb.score.MPV'] = self.posterior[:,1]
		self.data['fb.score.father'] = self.posterior[:,1]
		maxI = self.posterior.argmax( axis=1 )
		self.data['fb.prediction'] = self.labels[maxI]
	
	def _pathF( self ):
		pass

class DecodeAll( DecodeViterbi, DecodeForwardBackward ):
	''' Viterbi and foward-backward decoding
	'''
	def run( self ):
		''' main function accessible by outside classes '''
		# Viterbi
		self._initializeV()
		self._fillV()
		self._pathV()
		# FB
		self._initializeF()
		self._fillF()
		self._pathF()
		return self.data
