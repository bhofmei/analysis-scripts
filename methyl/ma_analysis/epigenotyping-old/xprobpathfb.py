import pandas as pd
import numpy as np

class ProbPathFB:
	
	def __init__( self, df, transitions, isUniform=False ):
		self.labels = np.array(['mother','MPV','father'])
		self.transitions = np.array(transitions)[0:3,0:3]
		self.data = df
		self.emissions = self._getEmissions( self.data )
		self.states, self.size = self.emissions.shape
		self.uniform = isUniform
		#print('transitions')
		#print(self.transitions)
		#print('emissions', self.data['sample'].iloc[0] )
		#print(self.emissions[:,:5])
	
	def run( self ):
		self._initialize()
		self._fill()
		self._path()
		return self.data
	
	def _getEmissions( self, inData ):
		inDataMelt = pd.melt(inData, id_vars=['sample','bin','prediction'])
		inDataPivotTmp = inDataMelt.pivot(index='variable', columns='bin' ,values='value')
		inDataPivot = inDataPivotTmp.reindex( self.labels )
		#print(inDataPivot.head)
		return np.exp(inDataPivot.values)
	
	def _initialize( self ):
		self.forward = np.zeros((self.states, self.size+1))
		self.forward[:,0] = (1.0/self.states if self.uniform else np.array([0.25, 0.5, 0.25] ) )
		self.backward = np.zeros((self.states,self.size+1))
		self.backward[:,-1] = 1.0
		self.posterior = np.zeros((self.size,self.states))
	
	def _fill( self ):
		# fill forward
		for i in range(self.size):
			fRow = np.matrix(self.forward[:,i])
			self.forward[:,i+1] = fRow * np.matrix( self.transitions ) * np.matrix( np.diag(self.emissions[:,i]) )
			# normalize
			self.forward[:,i+1] = self.forward[:,i+1] / np.sum(self.forward[:,i+1])
		#print('forward')
		#print( self.forward[:,:5])
		# fill backward
		for i in range(self.size, 0, -1):
			bCol = np.matrix( self.backward[:,i]).transpose()
			self.backward[:,i-1] = ( np.matrix(self.transitions) * np.matrix( np.diag(self.emissions[:,i-1])) * bCol ).transpose()
			# normalize
			self.backward[:,i-1] = self.backward[:,i-1] / np.sum( self.backward[:,i-1])
		#print('backward')
		#print( self.backward[:,-5:])
		# combine
		tmpPosterior = np.zeros(( self.states, self.size))
		tmpPosterior = np.array( self.forward[:,1:] ) * np.array( self.backward[:,:-1] )
		tmpPosterior = tmpPosterior / np.sum( tmpPosterior, 0)
		self.posterior = np.transpose(tmpPosterior)
		#print(self.posterior.shape)
	
	def _path( self ):
		self.data['score'] = self.posterior.max( axis=1 )
		maxI = self.posterior.argmax( axis=1 )
		self.data['iprediction'] = self.labels[maxI]
	

	'''def _fillForward( self ):
		# iterate through rows/bins
		for i in range(self.size):
			# iterate through states
			for j in range(3):
				emissionProb = self.data[[self.labels[j]]].iloc[i]
				self.forwardMat[i][j] = self._computeFwdScore( i, j, emissionProb )
			# end for j
			scaleTmp = sum(self.forwardMat[i])
			self.forwardMat[i] = [self.forwardMat[i][j] / scaleTmp for j in range(3) ]
		# end for i
	# normalize forward mat once we have filled the whole table
	
	def _computeFwdScore( self, i, j, prob ):
		tmpSum = 0
		uniform = [1.0/3,1.0/3,1.0/3]
		nonuniform = [0.25, 0.5, 0.25]
		for k in range(3):
			if i == 0:
				tmpSum += ( uniform[k] if self.uniform else nonuniform[k] ) * self.transMat[k][j]
			else:
				tmpSum += self.forwardMat[i-1][k] * self.transMat[k][j]
		return tmpSum * prob
	
	def _fillBackward( self ):
		# iterate reverse through rows/bins
		for i in range(self.size-1, -1, -1):
			# iterate through states
			for j in range(3):
				self.backwardMat[i][j] = self._computeBwdScore( i, j )
			# end for j
			# normalize
			scaleTmp = sum(self.backwardMat[i])
			self.backwardMat[i] = [self.backwardMat[i][j] / scaleTmp for j in range(3) ]
		# end for i
	
	def _computeBwdScore( self, i, j ):
		tmpSum = 0
		if i == self.size - 1:
			uniform = [1.0/3,1.0/3,1.0/3]
			nonuniform = [0.25, 0.5, 0.25]
			init  = ( uniform if self.uniform else nonuniform )
			for k in range(3):
				tmpSum += self.transMat[j][k] * init[k]
		else:
			for k in range(3):
				em = self.data[[self.labels[k]]].iloc[i+1]
				tr = self.transMat[j][k]
				tmpSum += self.backwardMat[i+1][k] * em * tr
		return tmpSum
	
	def _fillProbMat( self ):
		for i in range(self.size):
			for j in range(3):
				self.posteriorMat[i][j] = self.forwardMat[i][j] * self.backwardMat[i][j]
			scaleTmp = sum( self.posteriorMat[i] )
			self.posteriorMat[i] = [ self.posteriorMat[i][j] / scaleTmp for j in range(3) ]
			# end for j
		# end for i
	
	def getPath( self ):
		# add columns to data
		self.data['score'] = 0
		ss = np.nonzero( self.data.columns.values == 'score' )[0][0]
		self.data['iprediction'] = 'NA'
		si = np.nonzero( self.data.columns.values == 'iprediction' )[0][0]
		for i in range(self.size):
			vals = self.posteriorMat[i]
			maxS = max( vals )
			maxP = vals.index(maxS)
			self.data.iloc[i,ss] = maxS
			self.data.iloc[i,si] = self.labels[maxP]
	
	def run( self ):
		self._initializeMat()
		self._fillForward()
		self._fillBackward()
		self._fillProbMat()
		self.getPath()
		return self.data '''
		
