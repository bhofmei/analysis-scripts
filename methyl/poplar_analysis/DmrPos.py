class DmrPos:
	def __init__(self, chrm, start, end, label=''):
		self.chrm = chrm
		self.start = start
		self.end = end
		self.label = label
		self.posAr = []
		self.negAr = []
	
	def __str__( self ):
		outStr = '{:s}\t{:d}\t{:d}\t{:s}'.format( self.chrm, self.start, self.end, self.label)
		outStr += '\t' + ','.join([ str(x) for x in self.posAr])
		outStr += '\t' + ','.join([str(y) for y in self.negAr])
		return outStr
	
	def __eq__( self, other ):
		return self.chrm == other.chrm and self.start == other.start and self.end == other.end
	
	def __lt__( self, other ):
		if self.chrm == other.chrm:
			if self.start == other.start:
				return self.end < other.end
			return self.start < other.start
		return self.chrm < other.chrm
	
	def setPos( self, posAr, negAr ):
		self.posAr = posAr
		self.negAr = negAr
	
	def getCoord( self ):
		return self.chrm, self.start, self.end
	
	def getInfo ( self ):
		return self.chrm, self.start, self.end, self.label
