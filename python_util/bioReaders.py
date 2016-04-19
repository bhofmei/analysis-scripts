# class for reading different input file types
import sys, os

class BioReader:

	def __init__(self, fileStr, fileType=None):
	
		self.inFileStr = fileStr
		self.fileType = fileType

	def __repr__( self ):
		if fileType != None:
			return '{:s}: {:s}'.format( self.fileType, self.inFileStr )
		else:
			return self.inFileStr
	
	def exists( self ):
		return os.path.isfile( self.inFileStr )
	
	### GFF readers ###
	def readGFF_GenesDict( self ):
	
		gffFile = open (self.inFileStr, 'r' )
		gffDict = {}
	
		for line in gffFile:
			if line.startswith('#'):
				continue
			lineAr = line.rstrip().split('\t')
			# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
			# (6) strand (7) ? (8) notes
			start = int(lineAr[3])
			end = int(lineAr[4])
			chrm = lineAr[0]
			strand = lineAr[6]
		
			# only apply to type = gene
			if lineAr[2] == "gene":
				name = getGeneName( lineAr[8] )
				# put into geneDict
				gffDict[name] = (chrm, start, end, strand)
	
		gffFile.close()
		return gffDict

	def getGeneName(self, notesStr):
		search = "Name="
		index = notesStr.find(search)
		adIndex = index + len(search)
		endIndex = notesStr[adIndex:].find(';')
		if endIndex == -1:
			return notesStr[adIndex:]
		else:
			return  notesStr[adIndex:endIndex+adIndex]

	def readGFF_GenesAr( self ):
		gffFile = open (self.inFileStr, 'r' )
		gffAr = []
	
		for line in gffFile:
			if line.startswith('#'):
				continue
			lineAr = line.rstrip().split('\t')
			# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
			# (6) strand (7) ? (8) notes
			start = int(lineAr[3])
			end = int(lineAr[4])
			chrm = lineAr[0]
			strand = lineAr[6]
		
			# only apply to type = gene
			if lineAr[2] == "gene":
				gffAr += [(chrm, start, end, strand)]
	
		gffFile.close()
		return gffAr

	def readBedFR( self, bedFileStr ):
		'''
			creates a dictionary with each scaffold and those dictionaries are a 
			dictionary
			for the positions of coordinates with frequency count from the bed file
			{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
			return the dictionary
		'''
		bedFile = open( bedFileStr, 'r' )
		bedDict = {}
		count = 0
	
		# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
		for line in bedFile:
			lineAr =line.rstrip().split('\t')
			for pos in range( int(lineAr[1])+1, int(lineAr[2])+2):
				try:
					curScaf = lineAr[0]
					pos = int( lineAr[1] ) + 1
				except ValueError:
					pass
				else:
					count += 1
					# no dictionary for scaffold
					if bedDict.get(curScaf) == None:
						bedDict[curScaf] = {pos:1}
					# dictionary for scaffold but position not included
					elif bedDict.get(curScaf).get(pos) == None:
						bedDict[curScaf][pos] = 1
					# dictionary for scaffold and position included
					else:
						bedDict[curScaf][pos] += 1
	
		bedFile.close()
		return bedDict, count
	
