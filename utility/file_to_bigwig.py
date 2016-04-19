import sys, math, glob, multiprocessing, subprocess, os

# Usage: fileToBigWig.py [-keep] [-strand] [-scale] [-union] <bam_file | bed_file> <chrm_file>

def processInputs( bedFileStr, chrmFileStr, keepTmp, isStrand, isScale, isUnion ):
	print( 'Input file: {:s}\nChromosome sizes file: {:s}\nKeep temporary files: {:s}\nStrand-specific: {:s}\nScale by library size: {:s}'.format( bedFileStr, chrmFileStr, str( keepTmp ), str( isStrand ), str( isScale) ) )
	if isStrand:
		print( 'Combine strand-specific output: {:s}\n'.format( str(isUnion) ) )
		
	baseDir = os.path.dirname( bedFileStr )
	ind = bedFileStr.rfind('.')
	baseName = bedFileStr[:ind]
	rmFile = []
	
	# check for bam file -> convert to bed
	if bedFileStr.endswith( '.bam' ):
		print( 'Converting BAM to bed...' )
		bamFileStr = bedFileStr
		bedFileStr = '{:s}.bed'.format( baseName)
		command = 'bedtools bamtobed -i {:s} > {:s}'.format( bamFileStr, bedFileStr )
		subprocess.call( command, shell=True )
		rmFile += [ bedFileStr ]
		
	# get scale value if necessary
	scaleVal = 1
	if isScale:
		scaleVal = getScaleValue( bedFileStr )
	
	# bed to bedGraph
	if isUnion:
		bedGraphFileAr = convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal )
		rmFile += bedGraphFileAr
		bedGraphFile = uniteBedGraph( bedGraphFileAr[0], bedGraphFileAr[1], chrmFileStr, baseName )
		subAr = [ '_union' ]
	elif isStrand:
		bedGraphFile = convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal )
		subAr = [ '_plus', '_minus' ]
	else:
		bedGraphFile = convertToBedGraph (bedFileStr, chrmFileStr, baseName, scaleVal )
		subAr = [ '' ]
	
	rmFile += bedGraphFile
	
	# convert to bigwig
	bigWigFile = convertToBigWig( bedGraphFile, chrmFileStr, baseName, subAr )
	
	if keepTmp == False:
		print( 'Removing temporary files...' )
		for f in rmFile:
			os.remove(f)	
	
	print( '\nOutput written to {:s}'.format( ' & '.join(bigWigFile) ) )

def getScaleValue( bedFileStr ):
	command = 'wc -l {:s}'.format( bedFileStr )
	readCountStr = subprocess.check_output( command, shell=True, universal_newlines = True )
	return float( readCountStr.split()[0] ) / 1000000

def convertToBedGraph( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating bedGraph...' )
	bedGraphFile = '{:s}.bedGraph'.format( baseName )
	command = 'bedtools genomecov -bga -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFile )
	subprocess.call( command, shell=True )
	return [bedGraphFile]

def convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating strand-specific bedGraph...' )
	bedGraphFilePlus = '{:s}_scale_plus.bedGraph'.format( baseName )
	bedGraphFileMinus = '{:s}_scale_minus.bedGraph'.format( baseName )
	
	command = 'bedtools genomecov -bga -strand + -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFilePlus )
	subprocess.call( command, shell=True )
	command = 'bedtools genomecov -bga -strand - -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFileMinus+'.tmp' )
	subprocess.call( command, shell=True )
	# negative strand -> negative scores
	negativeBedGraph( bedGraphFileMinus, bedGraphFileMinus + '.tmp' )
	os.remove( bedGraphFileMinus + '.tmp' )
	return [bedGraphFilePlus, bedGraphFileMinus]

def negativeBedGraph( bedGraphFileMinus, bedGraphTmp ):
	
	tmpFile = open( bedGraphTmp, 'r' )
	bedGraphFile = open( bedGraphFileMinus, 'w' )
	
	for line in tmpFile:
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) > 3:
			lineAr[3] = '-'+lineAr[3]
		bedGraphFile.write( '{:s}\n'.format( '\t'.join(lineAr) ) )
	tmpFile.close()
	bedGraphFile.close()

def uniteBedGraph( bedGraphPlus, bedGraphMinus, chrmFileStr, baseName ):
	print( 'Uniting bedGraphs...' )
	bedGraphFile = '{:s}_stranded.bedGraph'.format( baseName )
	command = 'bedtools unionbedg -empty -g {:s} -i {:s} {:s} > {:s}'.format( chrmFileStr, bedGraphPlus, bedGraphMinus, bedGraphFile + '.tmp' )
	subprocess.call( command, shell=True )
	sumBedGraph( bedGraphFile, bedGraphFile + '.tmp' )
	os.remove( bedGraphFile + '.tmp' )
	return [ bedGraphFile ]
	
def sumBedGraph( bedGraphFile, bedGraphTmp ):
	
	bedGraph = open( bedGraphFile, 'w' )
	tmpFile = open( bedGraphTmp, 'r' )
	for line in tmpFile:
		lineAr = line.rstrip().split('\t' )
		# (0) chrm (1) start (2) end (3+) samples
		lNum = [ float(x) for x in lineAr[3:] ]
		bedGraph.write( '{:s}\t{:.2f}\n'.format( '\t'.join(lineAr[0:3]), sum(lNum) ) )
	tmpFile.close()
	bedGraph.close()
	
def convertToBigWig( bedGraphFileAr, chrmFileStr, baseName, subAr ):
	
	print( 'Converting to bigWig...' )
	bwFiles = [ '{:s}{:s}.bw'.format(baseName, x) for x in subAr ]
	for i in range(len(bedGraphFileAr)):
		command = 'bedGraphToBigWig {:s} {:s} {:s}'.format( bedGraphFileAr[i], chrmFileStr, bwFiles[i] )
		subprocess.call( command, shell=True )
	return bwFiles

def parseInputs( argv ):
	keepTmp = False
	isStrand = False
	isScale = False
	isUnion = False
	startInd = 0
	for i in range(min(4,len(argv)-2)):
		if argv[i] == '-keep':
			keepTmp = True
			startInd += 1
		elif argv[i] == '-strand':
			isStrand = True
			startInd += 1
		elif argv[i] == '-scale':
			isScale = True
			startInd += 1
		elif argv[i] == '-union':
			isUnion = True
			startInd += 1
	# can't have union and not strand
	if isUnion and isStrand == False:
		print( 'ERROR: union should only be used with strand' )
		exit()
	bedFileStr = argv[startInd]
	chrmFileStr = argv[startInd+1]
	
	processInputs( bedFileStr, chrmFileStr, keepTmp, isStrand, isScale, isUnion )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: fileToBigWig.py [-keep] [-strand] [-scale] [-union] <bam_file | bed_file> <chrm_file>")
	else:
		parseInputs( sys.argv[1:] )
