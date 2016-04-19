# python histone modification colors definitions

class hmColors:
	
	def __init__(self):
		self.colorDict = { 'h2az':'#ee7600', 'h3':'#8b7765', 'h3k4m3':'#9a32cd',
		'h3k9m2':'#228b22', 'h3k56ac':'#ee1289', 'input':'#708090',
		'h3k36m3':'#ee2c2c', 'h3k27m3':'#3a5fcd', 'h3t32':'#008b8b',
		'h2bs112':'#6d32c3', 'genes':'#daa520', 'gene':'#daa520', 
		'rnas':'#18A071', 'rna':'#18A071', 'tes':'#77158D',
		'transposons':'#77158D', 'mcg':'#b03060', 'mchg':'#2e8b57',
		'mchh':'#1e90ff', 'cg':'#b03060', 'chg':'#2e8b57',
		'chh':'#1e90ff','h3k27m3':'#617ed7','h3t32':'#32a2a2',
		'h3k36m1':'#AE2020','h3k36m2':'#D42727', 'h3k4m1':'#6A228D',
		'h3k4m2':'#872CB3','sdg7':'#2e8b57','basej':'#228b22',
		'white':'#ffffff' }
	
	def __repr__(self):
		return 'histone modification colors'
	
	def getColorStr( self, typeStr ):
		outStr = self.colorDict.get( typeStr.lower() )
		return outStr
	
	def getColorAr( self, typeStr ):
		outStr = self.colorDict.get( typeStr.lower() )
		if outStr == None:
			return [0,0,0]
		else:
			r,g,b = outStr[1:3], outStr[3:5], outStr[5:7]
			y = [ int(x,base=16) for x in [r,g,b] ]
			return y
	
	def getColorNormAr( self, typeStr ):
		outStr = self.colorDict.get( typeStr.lower() )
		if outStr == None:
			return [0,0,0]
		else:
			r,g,b = outStr[1:3], outStr[3:5], outStr[5:7]
			y = [ int(x,base=16) for x in [r,g,b] ]
			yy = [ float(x)/255 for x in y ]
			return yy
