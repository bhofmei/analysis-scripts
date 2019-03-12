# probability path that takes into account length of block

import pandas as pd
import numpy as np
from probpathsimple import ProbPathSimple

class ProbPathLength( ProbPathSimple ) :
	def _fillMat( self ):
		labels = ['mother', 'MPV', 'father' ]
		count = 0
		# loop through row/bin
		for i in range(self.size):
			for j in range(3):
				prob = self.data[labels[j]].iloc[i]
				maxScore, maxPos, count = self._computeScore( i, j, prob, count )
				self.probMat[i][j] = maxScore
				self.traceMat[i][j] = maxPos
			# end for j
		# end for i
	
	def _computeScore( self, i, j, prob, count ):
		# note self.transMat is actually a dictionary
		lbls = ['mother', 'MPV', 'father']
		scores = [prob]*3
		for k in range(3):
			# need to compute appropriate transitions
			if j == k:
				a = self.transMat[lbls[j]][lbls[j]]
				l = min( len(a)-1, count )
				trans = a[l]
			else:
				trans = self.transMat[lbls[j]][lbls[k]]
			scores[k] += trans + (0 if i == 0 else self.probMat[i-1][k])
			#print( i,j,k, prob, trans,(0 if i == 0 else self.probMat[i-1][k]), scores )
		print(scores)
		maxScore = max(scores)
		maxIndex = scores.index(maxScore)
		newCount = (count+1 if maxIndex == j else 0)
		maxPos = (None if i == 0 else maxIndex)
		return maxScore, maxPos, newCount
		
