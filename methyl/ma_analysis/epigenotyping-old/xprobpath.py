# Optimal path types
import numpy as np
import pandas as pd

class ProbPath:
	
	def __init__( self, df, transitions, isUniform ):
		self.labels = np.array(['mother', 'MPV', 'father'])
		self.data = df
		self.transitions = transitions
		self.emissions = self._getEmissions( df )
		self.states, self.size = self.emissions.shape
		#print(self.states,self.size)
		self.uniform = isUniform
	
	def _getEmissions( self, inData ):
		inDataMelt = pd.melt(inData, id_vars=['sample','bin','prediction'])
		inDataPivotTmp = inDataMelt.pivot(index='variable', columns='bin' ,values='value')
		inDataPivot = inDataPivotTmp.reindex( self.labels )
		return inDataPivot.values
	
class ProbPathViterbi( ProbPath ):

	def run( self ):
		# initialize
		self._initializeV()
		# fill matrix
		self._fillV()
		# get paths
		self._pathV()
		return self.data
	
	def _initializeV( self ):
		# get log of transitions
		self.log_transitions = np.log( self.transitions )
		self.probabilities = np.zeros( (self.size, self.states) )
		self.traceback = np.zeros( (self.size, self.states), dtype=np.int8 )
		nbin = self.size
		
	def _fillV( self ):
		'''for i in range(self.size):
			if i == 0:
				prev = np.log(np.array( ([1.0/3]*3 if self.uniform else [0.25,0.5,0.25]) ) )
			else:
				prev = self.probabilities[i-1,:]
			for j in range(self.states):
				em = self.emissions[j,i]
				scores = em + (prev + self.log_transitions[:,j])
				self.probabilities[i,j] = scores.max()
				self.traceback[i,j] = ( -1 if i == 0 else scores.argmax()  )
			# end for j
		# end for i'''
		# loop through rows/bin
		for i in range(self.size):
			# loop through columns/states
			for j in range(3):
				prob = self.emissions[j,i]
				maxS, maxP = self._computeScore( i, j, prob )
				#print( prob, maxS )
				self.probabilities[i,j] = maxS
				self.traceback[i,j] = maxP
			# end for j
		# end for i
	
	def _computeScore( self, i, j, prob ):
		scores = [prob]*3
		for k in range(3):
			#scores[k] *= ( 1 if i == 0 else self.probMat[i-1][k]) * self.transMat[k][j]
			# change this to use log scores
			scores[k] += ( 0 if i == 0 else self.probabilities[i-1,k]) + self.log_transitions[k,j]
		#print( scores )
		maxS = max( scores )
		maxI = scores.index( maxS )
		maxPos = ( -1 if i == 0 else maxI )
		return maxS, maxPos
	
	
	def _pathV( self ):
		# add columns to dataframe for output
		#self.data['score'] = 0
		self.data['vit.score.mother'] = self.probabilities[:,0]
		self.data['vit.score.MPV'] = self.probabilities[:,1]
		self.data['vit.score.father'] = self.probabilities[:,2]
		self.data['vit.prediction'] = 'NA'
		vals = self.probabilities[self.size-1]
		maxS = vals.max()
		nextJ = vals.argmax()
		for i in range(self.size-1, -1, -1):
			nextJ = self._tracebackHelper( i, nextJ )
			if nextJ == -1:
				break
		# end for i
	
	def _tracebackHelper( self, i, j ):
		#sj = np.nonzero(self.data.columns.values == 'score')[0][0]
		si = np.nonzero(self.data.columns.values == 'vit.prediction')[0][0]
		#score = self.probabilities[i,j]
		label = self.labels[j]
		#self.data.iloc[i,sj] = score
		self.data.iloc[i,si] = label
		return self.traceback[i,j]
	
class ProbPathFB( ProbPath ):
	
	def run( self ):
		# initialize
		self._initializeFB()
		# fill matrix
		self._fillFB()
		# get paths
		self._pathFB()
		return self.data

	def _initializeFB( self ):
		self.prob_emissions = np.exp( self.emissions )
		self.forward = np.zeros((self.states, self.size+1))
		self.forward[:,0] = (1.0/self.states if self.uniform else np.array([0.25, 0.5, 0.25] ) )
		self.backward = np.zeros((self.states,self.size+1))
		self.backward[:,-1] = 1.0
		self.posterior = np.zeros((self.size,self.states))
		
	def _fillFB( self ):
		# fill forward
		for i in range(self.size):
			fRow = np.matrix(self.forward[:,i])
			self.forward[:,i+1] = fRow * np.matrix( self.transitions ) * np.matrix( np.diag(self.prob_emissions[:,i]) )
			# normalize
			self.forward[:,i+1] = self.forward[:,i+1] / np.sum(self.forward[:,i+1])
			
		# fill backward
		for i in range(self.size, 0, -1):
			bCol = np.matrix( self.backward[:,i]).transpose()
			self.backward[:,i-1] = ( np.matrix(self.transitions) * np.matrix( np.diag(self.prob_emissions[:,i-1])) * bCol ).transpose()
			# normalize
			self.backward[:,i-1] = self.backward[:,i-1] / np.sum( self.backward[:,i-1])
			
		# combine
		tmpPosterior = np.zeros(( self.states, self.size))
		tmpPosterior = np.array( self.forward[:,1:] ) * np.array( self.backward[:,:-1] )
		tmpPosterior = tmpPosterior / np.sum( tmpPosterior, 0)
		self.posterior = np.transpose(tmpPosterior)
	
	def _pathFB( self ):
		self.data['fb.score.mother'] = self.posterior[:,0]
		self.data['fb.score.MPV'] = self.posterior[:,1]
		self.data['fb.score.father'] = self.posterior[:,2]
		maxI = self.posterior.argmax( axis=1 )
		self.data['fb.prediction'] = self.labels[maxI]

class ProbPathAll( ProbPathViterbi, ProbPathFB ):
	
	def run( self ):
		# Viterbi
		self._initializeV()
		self._fillV()
		self._pathV()
		# FB
		self._initializeFB()
		self._fillFB()
		self._pathFB()
		return self.data
