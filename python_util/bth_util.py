# useful functions for custom scripts
# bth_util.py

import sys, os

def binSizeToStr( binSize ):
	'''
		converts a bin size (bp) to a formatted string that is easier to read
		only converts when bin size divides reasonably easily
		i.e. 1,000,000 = 1mbp, 2,500 = 2.5kbp, 3,452 = 3452bp
	'''
	gb = 1000000000
	mb = 1000000
	kb = 1000
	
	# check gbp
	if ( binSize // gb ) > 0:
		x = binSize / gb
		if x % 1 == 0:
			return '{:d}gbp'.format( int(x) )
		elif x % 1 == 0.5:
			return '{:.1f}gbp'.format( x )
		return '{:d}bp'.format( binSize )
	# check mbp
	elif ( binSize // mb ) > 0:
		x = binSize / mb
		if x % 1 == 0:
			return '{:d}mbp'.format( int(x) )
		elif x % 1 == 0.5:
			return '{:.1f}mbp'.format( x )
		return '{:d}bp'.format( binSize )
	# check kbp
	elif ( binSize // kb ) > 0:
		x = binSize / kb
		if x % 1 == 0:
			return '{:d}kbp'.format( int(x) )
		elif x % 1 == 0.5:
			return '{:.1f}kbp'.format( x )
		return '{:d}bp'.format( binSize )
	# bp
	return '{:d}bp'.format( binSize ) 

def fileBaseName( fileStr ):
	bName = os.path.basename( fileStr )
	rInd = bName.rfind( '.' )
	if rInd != -1:
		return bName[:rInd]
	return bName
	
def rgb2hex( r, g, b ):
	sR = hex(r)[2:]
	sG = hex(g)[2:]
	sB = hex(b)[2:]
	return sR + sG + sB

def hex2rgb( hStr ):
	hStr = hStr.replace('#', '' )
	if len(hStr) > 6:
		print( 'ERROR: {:s} is not valid -- too long' )
		exit()
	hR = hStr[0:2]
	hG = hStr[2:4]
	hB = hStr[4:6]
	iR = int( hR, base = 16 )
	iG = int( hG, base = 16 )
	iB = int( hB, base = 16 )
	return ( iR, iG, iB )
	
def strToDistance( inStr ):
	
	inStr = inStr.upper()
	# kilo basepair
	if inStr.endswith( 'K' ):
		return int( inStr[:-1] ) * 1000
	elif inStr.endswith( 'KB' ):
		return int( inStr[:-2] ) * 1000
	elif inStr.endswith( 'KBP' ):
		return int( inStr[:-3] ) * 1000
	# mega basepair
	elif inStr.endswith( 'M' ):
		return int( inStr[:-1] ) * 1000000
	elif inStr.endswith( 'MB' ):
		return int( inStr[:-2] ) * 1000000
	elif inStr.endswith( 'MBP' ):
		return int( inStr[:-3] ) * 1000000
	# giga basepair
	elif inStr.endswith( 'G' ):
		return int( inStr[:-1] ) * 1000000000
	elif inStr.endswith( 'GB' ):
		return int( inStr[:-2] ) * 1000000000
	elif inStr.endswith( 'GBP' ):
		return int( inStr[:-3] ) * 1000000000
	# basepair
	elif inStr.endswith( 'BP' ):
		return int( inStr[:-2] )
	else:
		try:
			return int( inStr )
		except ValueError:
			return False
