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
		elif (x*10) % 1 == 0.0:
			return '{:.1f}gbp'.format( x )
		elif (x*100) % 1 == 0.0:
			return '{:.2f}gbp'.format( x )
		return '{:d}bp'.format( binSize )
	# check mbp
	elif ( binSize // mb ) > 0:
		x = binSize / mb
		if x % 1 == 0:
			return '{:d}mbp'.format( int(x) )
		elif (x*10) % 1 == 0.0:
			return '{:.1f}mbp'.format( x )
		elif (x*100) % 1 == 0.0:
			return '{:.2f}mbp'.format( x )
		return '{:d}bp'.format( binSize )
	# check kbp
	elif ( binSize // kb ) > 0:
		x = binSize / kb
		if x % 1 == 0:
			return '{:d}kbp'.format( int(x) )
		elif (x*10) % 1 == 0.0:
			return '{:.1f}kbp'.format( x )
		elif (x*100) % 1 == 0.0:
			return '{:.2f}kbp'.format( x )
		return '{:d}bp'.format( binSize )
	# bp
	return '{:d}bp'.format( binSize ) 

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

def fileBaseName( fileStr ):
	bName = os.path.basename( fileStr )
	rInd = bName.rfind( '.' )
	if rInd != -1:
		return bName[:rInd]
	return bName

''' color conversion functions '''	
def rgb2hex( r, g, b, inS=False ):
	sR = hex(r)[2:]
	if len(sR) == 1:
		sR = '0'+sR
	sG = hex(g)[2:]
	if len(sG)==1:
		sG = '0' + sG
	sB = hex(b)[2:]
	if len(sB) == 1:
		sB = '0' + sB
	return ('#' if inS else '') + sR + sG + sB

def rgb2dec( r, g, b ):
	fR = float( r / 255 )
	fG = float( g / 255 )
	fB = float( b / 255 )
	return ( fR, fG, fB )

def rgb2decf( r, g, b ):
	iR, iG, iB = rgb2dec (r, g, b)
	return decFormat( iR, iG, iB )

def decFormat( r, g, b ):
	return '({:.3f}, {:.3f}, {:.3f})'.format( r, g, b )

def dec2hex( r, g, b ):
	iR, iG, iB = dec2rgb( r, g, b )
	return rgb2hex( iR, iG, iB )

def dec2rgb( r, g, b ):
	iR = int( r*255 )
	iG = int( g*255 )
	iB = int( b*255 )
	return ( iR, iG, iB )

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

def hex2dec( hStr ):
	r, g, b = hex2rgb( hStr )
	return rgb2dec( r, g, b )

def hex2decf( hStr ):
	r, g, b = hex2rgb( hStr )
	return rgb2decf( r, g, b )

def rgb2hsl( r, g, b ):
	r, g, b = rgb2dec( r, g, b )
	
	cMax = max( r, g, b )
	iMax = [r, g, b].index( cMax )
	cMin = min( r, g, b )
	d = cMax - cMin
	l = (cMax+cMin) / 2
	
	if d == 0:
		h,s = 0,0
	else:
		s = d / ( 2.0-cMax-cMin  if l > 0.5 else cMax+cMin)
		if iMax == 0:
			h = (g-b)/d + (6 if g < b else 0)
		elif iMax == 1:
			h = (b-r)/d + 2
		elif iMax == 2:
			h = (r-g) / d + 4
		h /= 6
	return (h,s,l)
	
def hex2hsl( hStr ):
	r, g, b = hex2rgb( hStr )
	return rgb2hsl( r, g, b )

def hsl2dec( h, s, l ):
	# h should be angle, convert to dec
	h /= 360.0
	if s == 0:
		dR, dG, dB = [l]*3
	else:
		q = l*(1+s) if l < 0.5 else l + s - l*s
		p = 2*l - q
		dR = hue2rgb( p, q, h+1.0/3 )
		dG = hue2rgb( p, q, h )
		dB = hue2rgb( p, q, h-1.0/3 )
	return ( dR, dG, dB )

def hsl2rgb( h, s, l ):
	fR, fG, fB = hsl2dec( h, s, l )
	return dec2rgb(fR, fG, fB)

def hsl2hex( h, s, l ):
	# h should be angle
	r, g, b = hsl2rgb( h, s, l )
	return rgb2hex( r, g, b, True )
	
def hue2rgb(p, q, t):
	if t < 0:
		t += 1
	if t > 1:
		t -=1
	if t < 1.0/6:
		return p + (q-p) * 6 * t
	elif t < 0.5:
		return q
	elif t < 2.0/3.0:
		return p + (q-p) * (2.0/3.0 - t) * 6
	else:
		return p
	
def hexSplitComplimentary( hStr ):
	# convert to hsl
	h, s, l = hex2hsl( hStr )
	h = h *360.0
	nHs = [ h+x for x in [-150, 0, 150] ]
	nHs = [ (360+x if x < 0 else x) for x in nHs ]
	nVals = [ (nH, s, l) for nH in nHs ]
	outAr = []
	for i in range(len(nVals)):
		x,y,z = nVals[i]
		hexx = hsl2hex( x,y,z )
		outAr += [ hexx ]
	return outAr

def hex2dec_str( hStr ):
	r, g, b = hex2dec( hStr )
	return '({:.4f}, {:.4f}, {:.4f})'.format( r, g, b )
