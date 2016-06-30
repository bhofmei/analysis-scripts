# different classes for different types of biological data files
import sys, os, pickle, math, gzip

class FileBio:
	
	def __init__( self, inFileStr ):
		self.fileStr = inFileStr
		self.isZip = self.__isGZip( inFileStr )
	
	def __str__( self ):
		return os.path.basename( self.fileStr )
	
	def __isGZip( self, fileStr ):
		if fileStr.endswith( '.gz' ) or fileStr.endswith( '.gzip' ):
			return True
		return False
	
	def fbOpen( self ):
		if self.isZip:
			return gzip.open( self.fileStr, 'rt' )
		else:
			return open( self.fileStr, 'r' )
	
	def fbBasename( self ):
		l = self.fileStr.replace('.gz', '' ).replace( '.gzip', '' )
		rInd = l.rfind( '.' )
		if rInd != -1:
			return l[:rInd]
		return l

##### GFF
class FileGFF( FileBio ):
	
	'''def __init__( self, inFileStr ):
		self.gffFileStr = inFileStr
	
	def __str__( self ):
		return os.path.basename( self.gffFileStr )'''

	def getGeneName(self, notesStr, rmPeriod=False):
		search = "Name="
		index = notesStr.find(search)
		adIndex = index + len(search)
		endIndex = notesStr[adIndex:].find(';')
		n = notesStr[adIndex:]
		if endIndex != -1:
			n = notesStr[adIndex:endIndex+adIndex]
		
		pInd = n.rfind( '.' )
		if rmPeriod and pInd != 1:
			return n[:pInd]
		else:
			return n

	def getGeneArray( self, chrm=False, te=False, subsetAr=None ):
		''' return list of genes
			genes stored as tuples (chrm, start, end, strand)
		'''
		gffFile = self.fbOpen( )
		gffAr = []
		returnChrm = chrm
		chrmFormat = None
	
		for line in gffFile:
			if line.startswith('#'):
				continue
			lineAr = line.rstrip().split('\t')
			# (0) scaffold (1) organism (2) type (3) start (4) end (5) ? 
			# (6) strand (7) ? (8) notes
			chrm = lineAr[0]
			start = int(lineAr[3])
			end = int(lineAr[4])
			strand = lineAr[6]
			if chrmFormat == None:
				if chrm.isdigit():
					chrmFormat = 'digit'
				elif chrm.startswith( 'Chr0' ):
					chrmFormat = 'zeroed'
				elif chrm.startswith( 'Chr' ):
					chrmFormat = 'chr'
				else:
					chrmFormat = 'asis'
			
			if te and lineAr[2] in ['transposable_element', 'transposable_element_gene']:
				gffAr+= [(chrm, start, end, strand)]
			elif te == False and lineAr[2] == "gene":
				gName = self.getGeneName( lineAr[8] )
				if subsetAr == None or gName in subsetAr:
					gffAr+= [(chrm, start, end, strand)]
		gffFile.close()
		
		if returnChrm:
			return gffAr, chrmFormat
		else:
			return gffAr
	
	def getGeneDict( self, chrm=False, rmPeriod=False, includeTE=False ):
		''' return dictionary of genes by gene name
			value of dictionary is tuple (chrm, start, end, strand )
		'''
		gffFile = self.fbOpen()
		gffDict = {}
		returnChrm = chrm
		chrmFormat = None
		rP = rmPeriod
		
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
			
			if chrmFormat == None:
				if chrm.isdigit():
					chrmFormat = 'digit'
				elif chrm.startswith( 'Chr0' ):
					chrmFormat = 'zeroed'
				elif chrm.startswith( 'Chr' ):
					chrmFormat = 'chr'
				else:
					chrmFormat = 'asis'
			
			# only apply to type = gene
			if (includeTE and lineAr[2] in [ 'gene', 'transposable_element_gene', 'pseudogene' ]) or (lineAr[2] == 'gene' ):
				name = self.getGeneName( lineAr[8], rmPeriod = rP)
				# put into geneDict
				gffDict[name] = (chrm, start, end, strand)
	
		gffFile.close()
		if returnChrm:
			return gffDict, chrmFormat
		else:
			return gffDict

class FileBED( FileBio ):
	'''
		stores information about the start or midpoint of read
	'''
	
	def getBedDict( self, middle=False, stranded=False ):
		''' return dictionary where key is chromosome, value is dictionary
			that dictionary has positions as key and value is number of reads
			that cover the position
			dictionary contains key '##counts##' that stores total number of
			reads in library
		'''
		# check for pickle file
		
		if middle and stranded:
			#print( '~~is middle')
			pickleFileStr = self.fbBasename() + '_mid_strand.pick'
		elif middle:
			pickleFileStr = self.fbBasename() + '_mid.pick'
		elif stranded:
			pickleFileStr = self.fbBasename() + '_strand.pick'
		else:
			pickleFileStr = self.fbBasename() + '.pick'
		
		fileExists = os.path.isfile( pickleFileStr )
		# if exists -> check modification time
		if fileExists:
			mTime = os.path.getmtime( self.fileStr )
			tTime = os.path.getmtime( pickleFileStr )
			fSize = os.path.getsize( pickleFileStr )
			if tTime < mTime or fSize < 100:
				print( '~updating {:s}'.format(os.path.basename(pickleFileStr)) )
				bedDict = self.__createPickleFile( pickleFileStr, middle )
			else:
				print( '~loading {:s}'.format(os.path.basename(pickleFileStr)) )
				bedDict = pickle.load( open(pickleFileStr, 'rb') )
		# does not exist -> create
		else:
			bedDict = self.__createPickleFile( pickleFileStr, middle, stranded )
		# get read count
		counts = bedDict.get( '##counts##' )
		del bedDict['##counts##']
		return bedDict, counts

	def __createPickleFile( self, pickleFileStr, middle, stranded ):
		''' create the bed dictionary and save it to a pickle file
			return bedDict
		'''
		# create dictionary
		print( '~reading {:s}'.format( str(self) ) )
		bedDict, readCounts = self.__readBed( middle, stranded )
		# add count to dictionary as '##counts##'
		bedDict[ '##counts##' ] = readCounts
		# pickle
		print( '~creating {:s}'.format(os.path.basename(pickleFileStr)) )
		pickle.dump( bedDict, open( pickleFileStr, 'wb' ), protocol=3 )
		return bedDict
	
	def __readBed( self, middle, stranded ):
		'''
			creates a dictionary with each scaffold and those dictionaries are a 
			dictionary for the positions of coordinates with frequency count 
			from the bed file
			{scaf1:{pos1:1,pos2:4,pos3:1},scaf2{pos1:2}}
			return the dictionary
		'''
		#bedFile = open( self.bedFileStr, 'r' )
		bedFile = self.fbOpen()
		bedDict = {}
		readCount = 0
		
		# (0) scaffold (1) start (2) end (3) name (4) score (5) strand
		for line in bedFile:
			lineAr = line.rstrip().split('\t')
			try:
				curChrm = lineAr[0]
				end = int(lineAr[2])
				start = int(lineAr[1])
				if middle:
					pos = int( math.floor( (end - start + 1)/2.0 + start+1 ) )
				else:
					pos = int( lineAr[1] ) + 1
				if stranded:
					curChrm += lineAr[5]
			except ValueError:
				pass
			else:
				# no dictionary for chrm
				if bedDict.get(curChrm) == None:
					bedDict[curChrm] = {pos:1}
				# dictionary for chrm but position not included
				elif bedDict.get(curChrm).get(pos) == None:
					bedDict[curChrm][pos] = 1
				# dictionary for chrm and position included
				else:
					bedDict[curChrm][pos] += 1
				readCount += 1
	
		bedFile.close()
		return bedDict, readCount

class FileBED_FR( FileBio ):
	'''
		full read
		stores count information across the entire read
	'''
	
	def getBedDict( self ):
		''' return dictionary where key is chromosome, value is dictionary
			that dictionary has positions as key and value is number of reads
			that cover the position
			dictionary contains key '##counts##' that stores total bp in library
		'''
		# check for pickle file
		pickleFileStr = self.fbBasename() + '_fr.pick'
		
		fileExists = os.path.isfile( pickleFileStr )
		# if exists -> check modification time
		if fileExists:
			mTime = os.path.getmtime( self.fileStr )
			tTime = os.path.getmtime( pickleFileStr )
			fSize = os.path.getsize( pickleFileStr )
			if tTime < mTime or fSize < 100:
				print( '~updating {:s}'.format(os.path.basename(pickleFileStr)) )
				bedDict = self.__createPickleFile( pickleFileStr )
			else:
				print( '~loading {:s}'.format(os.path.basename(pickleFileStr)) )
				bedDict = pickle.load( open(pickleFileStr, 'rb') )
		# does not exist -> create
		else:
			bedDict = self.__createPickleFile( pickleFileStr )
		# get bp count
		counts = bedDict.get( '##counts##' )
		return bedDict, counts
	
	def __createPickleFile( self, pickleFileStr ):
		''' create the bed dictionary and save it to a pickle file
			return bedDict
		'''
		# create dictionary
		print( '~reading {:s}'.format( str(self) ) )
		bedDict, bpCounts = self.__readBed( )
		# add count to dictionary as '##counts##'
		bedDict[ '##counts##' ] = bpCounts
		# pickle
		print( '~creating {:s}'.format(os.path.basename(pickleFileStr)) )
		pickle.dump( bedDict, open( pickleFileStr, 'wb' ), protocol=4 )
		return bedDict
	
	def __readBed( self ):
		'''
			creates a dictionary with each scaffold and those dictionaries are a 
			dictionary for the positions of coordinates with frequency count from 
			the bed file
			i.e.: {chrm1:{pos1:1,pos2:4,pos3:1}, chrm2{pos1:2}}
			return the dictionary and total number of bp in library
		'''
		#bedFile = open( self.bedFileStr, 'r' )
		bedFile = self.fbOpen()
		bedDict = {}
		count = 0
	
		# (0) chrm (1) start (2) end (3) name (4) score (5) strand
		for line in bedFile:
			lineAr =line.rstrip().split('\t')
			for pos in range( int(lineAr[1])+1, int(lineAr[2])+2):
				try:
					curChrm = lineAr[0]
				except ValueError:
					pass
				else:
					count += 1
					# no dictionary for chrm
					if bedDict.get(curChrm) == None:
						bedDict[curChrm] = {pos:1}
					# dictionary for chrm but position not included
					elif bedDict.get(curChrm).get(pos) == None:
						bedDict[curChrm][pos] = 1
					# dictionary for chrm and position included
					else:
						bedDict[curChrm][pos] += 1
			# end for pos
	
		bedFile.close()
		return bedDict, count

class FileFPKM( FileBio ):
	'''
		reading files with FPKM values for genes
	'''
	'''def __init__( self, inFileStr ):
		self.fpkmFileStr = inFileStr
	
	def __str__( self ):
		return os.path.basename( self.fpkmFileStr )'''
	
	def getFPKMArray( self ):
		'''
			returns the array or genes ordered by fpkm
		'''
		# check if exists
		self.__checkFPKM( )
	
		fpkmFile = open( self.fileStr + ".txt", 'r' )
		fpkmAr = []
	
		for line in fpkmFile:
			line = line.rstrip()
			lineAr = line.split('\t')
			fpkmAr += [lineAr[0]]
		fpkmFile.close()
		return fpkmAr
	
	def getFPKMValueArray( self, subsetDict=None ):
		''' returns array of genes ordered by fpkm and array of the fpkms
			optional subsetDict returns only genes listed in that array 
			listed as the associated dictionary value
		'''
		# check if exists
		self.__checkFPKM()
		fpkmFile = open( self.fileStr + ".txt", 'r' )
		fpkmAr = []
		fpkmValueAr = []
		subsetAr = ( None if subsetDict == None else list(subsetDict.keys()) )
	
		for line in fpkmFile:
			line = line.rstrip()
			lineAr = line.split('\t')
			i = lineAr[0].rfind( '.' )
			if i != -1:
				lineAr[0] = lineAr[0][:i]
			# get values if needed
			if subsetAr == None and subsetDict == None:
				fpkmAr += [ lineAr[0] ]
				fpkmValueAr += [ lineAr[1] ]
			elif lineAr[0] in subsetAr:
				fpkmAr += [ subsetDict[lineAr[0]] ]
				fpkmValueAr += [ lineAr[1] ]
			fpkmAr += [lineAr[0]]
		fpkmFile.close()
		return fpkmAr

	def __checkFPKM( self ):
		'''
			checks of the self-created tab-deliminated fpkm file exists
			it it exists but is outdated, it is recreated
			if it doesn't exist, it is created
		'''
		newFpkmFile = self.fileStr + ".txt"
		fileExists = os.path.isfile( newFpkmFile )
	
		# if exists - check date modified and file size
		if fileExists:
			mTime = os.path.getmtime( self.fileStr )
			tTime = os.path.getmtime( newFpkmFile)
			fSize = os.path.getsize( newFpkmFile )
			if tTime < mTime or fSize < 100:
				print( '~updating {:s}'.format(newFpkmFile) )
				self.__createFPKMTextFile( newFpkmFile )
		# doesn't exist - create
		else:
			self.__createFPKMTextFile( newFpkmFile )
			print( '~creating {:s}'.format(newFpkmFile) )

	def __createFPKMTextFile( self, newFpkmFile ):
		'''
			creates the tab-deliminated FPKM file
			useful so gene order is the same across samples
		'''
		# fpkmDict set up as {fpkm:[gene1,gene2,...], fpkm2:[gene4],...}
		fpkmOutFile = open( newFpkmFile, 'w' )
		fpkmDict = self.__readFPKM( )
	
		for key in sorted(fpkmDict.keys(), reverse=True):
			genes = fpkmDict[key]
			# for each gene
			for gene in genes:
				fpkmOutFile.write( "{:s}\t{:f}\n".format( gene, key ) )
			
		fpkmOutFile.close()

	def __readFPKM( self ):
		'''
			creates a dictionary for genes by fpkm value
			fpkm value is key which points to an array of gene names
			return fpkm dictionary
		'''
		#fpkmFile = open( self.fpkmFileStr, 'r' )
		fpkmFile = self.fbOpen()
		fpkmDict = {}
	
		for line in fpkmFile:
			line = line.rstrip()
			lineAr = line.split('\t')
			# Adam's output
			# (0) locus (1) coverage (2) FPKM
			#print( len( lineAr ) )
			if len( lineAr ) == 3:
				# header
				if line.startswith( 'locus' ):
					continue
				fpkm = float( lineAr[2] )
				name = lineAr[0]
				if fpkmDict.get( fpkm ) == None:
					fpkmDict[fpkm] = [name]
				# case 2: fpkm in dict
				else:
					fpkmDict[fpkm] += [name]
			# Cufflinks output
			# (0) tracking_id (1) class_code (2) nearest_ref_id (3) gene_id 
			# (4) gene_short_name (5) tss_id (6) locus (7) length (8) coverage 
			# (9) FPKM (10) FPKM_conf_low (11) FPKM_conf_high (12) FPKM_status
			else:
				if lineAr[12] == 'OK':
					fpkm = float( lineAr[9] )
					name = lineAr[4]
					if name == '-':
						continue
					nameAr = name.split(',')
					for n in nameAr:
						# case 1: fpkm not in fpkmDict
						if fpkmDict.get( fpkm ) == None:
							fpkmDict[fpkm] = [n]
						# case 2: fpkm in dict
						else:
							fpkmDict[fpkm] += [n]
		return fpkmDict

class FileAllC_full( FileBio ):
	'''
		use when all chromosomes from allc are in the same file
		reading allC files and storing as a dictionary
	'''
	
	def getAllCDict( self, mtypes = None ):
		# check for pickle file
		pickleFileStr = self.fbBasename() + '.pick'
		
		fileExists = os.path.isfile( pickleFileStr )
		# if exists -> check modification time
		if fileExists:
			mTime = os.path.getmtime( self.fileStr )
			tTime = os.path.getmtime( pickleFileStr )
			fSize = os.path.getsize( pickleFileStr )
			if tTime < mTime or fSize < 100:
				print( '~updating {:s}'.format(os.path.basename(pickleFileStr)) )
				allcDict = self.__createPickleFile( pickleFileStr )
			else:
				print( '~loading {:s}'.format(os.path.basename(pickleFileStr)) )
				allcDict = pickle.load( open(pickleFileStr, 'rb') )
		# does not exist -> create
		else:
			allcDict = self.__createPickleFile( pickleFileStr )
		if mtypes == None:
			return allcDict
		if len(mtypes)==1:
			return allcDict[mtypes]
		s = sum([x in mtypes for x in allcDict.keys()])
		if s > 0:
			newDict = {}
			for x in mtypes:
				newDict[x] = allcDict[x]
			return newDict
		else:
			print( 'warning: allc dict empty; re-specify mtypes' )
			return False
	
	def __createPickleFile( self, pickleFileStr ):
		# create dictionary
		print( '~reading {:s}'.format( os.path.basename(self.fileStr) ) )
		allcDict = self.__readAllc( )
		# pickle
		print( '~creating {:s}'.format(os.path.basename(pickleFileStr)) )
		pickle.dump( allcDict, open( pickleFileStr, 'wb' ), protocol=4 )
		return allcDict
	
	def __readAllc( self, chrmFormat=None ):
		mTypes = [ 'CG', 'CHG', 'CHH', 'C' ]
		
		#allCFile = open( self.allcFileStr, 'r' )
		allCFile = self.fbOpen()
		allCDict = {}
	
		for m in mTypes:
			allCDict[m] = {}
	
		for line in allCFile:
			if line.startswith( 'chr\t' ):
				continue
			lineAr = line.rstrip().split('\t')
			# (0) chr (1) pos (2) strand (3) mc class (4) mc_count (5) total
			# (6) methylated
			if len( lineAr ) < 7 or lineAr[6].isdigit() == False:
				continue
			chrm = lineAr[0]
			# chrm format
			if chrmFormat == 'digit' and chrm.isdigit() == False:
				chrm = chrm.replace( 'Chr', '' )
			elif chrmFormat == 'zeroed' and chrm.isdigit():
				chrm = 'Chr{:02d}'.format( int(chrm) )
			elif chrmFormat == 'chr' and chrm.isdigit():
				chrm = 'Chr{:d}'.format( int( chrm ) )
			
			mLineType = self.findMethylType( lineAr[3] )
			if mLineType in mTypes:
				if allCDict[mLineType].get( chrm ) == None:
					allCDict[mLineType][chrm] = {}
				allCDict[mLineType][chrm][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
				
			if allCDict['C'].get( chrm ) == None:
				allCDict['C'][chrm] = {}
			allCDict['C'][chrm][int(lineAr[1])] = ( int(lineAr[4]), int( lineAr[5]) )
		# end for line
		allCFile.close()
		return allCDict
	
	def findMethylType( self, mc ):
		if mc.startswith( 'CG' ):
			return 'CG'
		elif mc.endswith( 'G' ):
			return 'CHG'
		elif mc == 'CNN':
			return 'CNN'
		else:
			return 'CHH'

class FileFASTAIndex( FileBio ):
	'''
		use for fasta index files, that is, files that end in .fa.fai
		this is a tab-deliminated file that includes chromosome names and
		lengths of the cooresponding fasta file
		columns: (0) chrm (1) length ...
	'''
	
	def getNames( self ):
		'''
			return array of chromosome names from file
		'''
		inFile = self.fbOpen()
		outAr = []
		for line in inFile:
			lineAr = line.rstrip().split( '\t' )
			outAr += [ lineAr[0] ]
		inFile.close()
		return outAr
	
	def getLengths( self, subset=None ):
		inFile = self.fbOpen()
		outDict = {}
		
		for line in inFile:
			lineAr = line.rstrip().split( '\t' )
			chrm = lineAr[0]
			if subset == None or chrm in subset:
				outDict[ chrm ] = int( lineAr[1] )
		inFile.close()
		return outDict

class FileFASTA( FileBio ):
	
	def getFastaDict( self ):
		'''
			return dictionary where key is chromosome, value is array
			array has FASTA sequence and is 1-based
		'''
		# check for pickle file
		pickleFileStr = self.fbBasename() + '.pick'
		
		fileExists = os.path.isfile( pickleFileStr )
		# if exists -> check modification time
		if fileExists:
			mTime = os.path.getmtime( self.fileStr )
			tTime = os.path.getmtime( pickleFileStr )
			fSize = os.path.getsize( pickleFileStr )
			if tTime < mTime or fSize < 100:
				print( '~updating {:s}'.format(os.path.basename(pickleFileStr)) )
				fastaDict = self.__createPickleFile( pickleFileStr )
			else:
				print( '~loading {:s}'.format(os.path.basename(pickleFileStr)) )
				fastaDict = pickle.load( open(pickleFileStr, 'rb') )
		# does not exist -> create
		else:
			fastaDict = self.__createPickleFile( pickleFileStr )
		return fastaDict
	
	def __createPickleFile( self, pickleFileStr ):
		''' create the fasta dictionary and save it to a pickle file
			return fastaDict
		'''
		# create dictionary
		print( '~reading {:s}'.format( str(self) ) )
		fastaDict = self.readFasta( )
		# pickle
		print( '~creating {:s}'.format(os.path.basename(pickleFileStr)) )
		pickle.dump( fastaDict, open( pickleFileStr, 'wb' ), protocol=3 )
		return fastaDict
	
	def readFasta( self ):
		''' 
			returns the created dictionary
		'''
		fastaFile = self.fbOpen()
		fastaDict = {}
		chrm = None
		seq = [' ']
		for line in fastaFile:
			line = line.rstrip()
			# chromosome headers
			if line.startswith('>'):
				# need to write old sequence
				if seq != [' '] and chrm != None:
					fastaDict[chrm] = seq
					seq = [' ']
				lineAr = line.split(' ')
				chrm = lineAr[0].replace('>', '')
		
			# sequence
			else:
				seqL = list(line.upper())
				seq += seqL
			
		# handle last chrm read
		if chrm not in fastaDict:
			fastaDict[chrm] = seq
		fastaFile.close()
		return fastaDict
		
