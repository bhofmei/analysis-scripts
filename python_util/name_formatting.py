# chromosome naming formatting

def getFormattingScheme( ):
	out =  'General (chrm, scaf, contig)\n'
	out += '\tu\tunderscore\n' 
	out += '\tz\tzeros prefixing 1\n' 
	out += '\te\tempty; just the digit (do not use with other options)\n' 
	out += 'Chromosomes:\tChr1 [default]\n'
	out += '\tl\tlowercase\n'
	out += '\th\tlong; \'chromosome\'\n' 
	out += 'Scaffolds:\tscaffold1 [default]\n' 
	out += '\tc\tcapitilize\n' 
	out += '\ts\tshort; \'scaf\'\n' 
	out += 'Contigs:\tcontig1 [default]\n' 
	out += '\tc\tcapitilize\n' 
	out += 'CLM:\n'
	out += '  Chloroplast:\tnot included [default]\n' 
	out += '\tC/c\tChrC/chrC\n' 
	out += '\tH/h\tChloroplast/chloroplast\n' 
	out += '  Lambda:\tnot included [default]\n' 
	out += '\tL/l\tChrL/chrL\n' 
	out += '\tA/a\tLambda/lambda\n' 
	out += '  Mitochondria:\tnot included [default]\n' 
	out += '\tM/m\tChrM/chrM\n' 
	out += '\tT/t\tChrMt/chrMt\n' 
	out += '\tI/i\tMitochondria/mitochondria\n'
	out += 'Other:\t\tas-is [default]\n'
	out += '\tc\tcapitalize\n'
	out += '\tl\tlowercase\n'
	out += '\th\tadd \'Chr\'\n'
	return out

def determineType( name ):
	
	n = name.lower()
	if n in ['chrc', 'chrm', 'chrl', 'chloroplast', 'mitochrondria', 'lambda', 'chrmt', 'c', 'm', 'l', 'mt' ]:
		return 'clm'
	elif n.startswith( 'chr' ):
		return 'chr'
	elif n.startswith( 'scaf' ):
		return 'scaf' 
	elif n.startswith( 'cont' ):
		return 'contig'
	elif n.isdigit():
		return 'chr'
	else:
		return 'other'

def checkEmpty( options, label ):
	isEmpty = 'e' in options.lower()
	if isEmpty and len(options) > 1:
		print( 'WARNING: do not specify empty with other options for {:s}...ignoring empty'.format(label) )
		options = options.lower().replace('e','')
	return options

def decodeChrmOptions( options ):
	opts = list(options.lower())
	isCap = 'l' not in opts
	isLong = 'h' in opts
	isUnSc = 'u' in opts
	isZero = opts.count( 'z' )
	isEmpty = 'e' in opts
	return isCap, isLong, isUnSc, isZero, isEmpty

def formatChrm( name, isCap, isLong, isUnSc, isZero, isEmpty ):
	out = ( 'Chr' if isCap else 'chr' )
	out += ('omosome' if isLong else '' )
	out += ( '_' if isUnSc else '' )
	if isEmpty:
		out = ''
	nname = name.lower().replace('chromosome','').replace('chr','').replace('_','')
	try:
		iname = int( nname )
		if isZero:
			out += '{number:0{width}d}'.format(width=(isZero+1), number=iname)
		else:
			out += str(iname)
	except ValueError:
		print( 'chrm error number (', nname, ')' )
		out += nname
	return out

def decodeScafOptions( options ):
	opts = list(options.lower())
	isCap = 'c' in opts
	isShort = 's' in opts
	isUnSc = 'u' in opts
	isZero = opts.count( 'z' )
	isEmpty = 'e' in opts
	return isCap, isShort, isUnSc, isZero, isEmpty

def formatScaf( name, isCap, isShort, isUnSc, isZero, isEmpty ):
	out = ('Scaf' if isCap else 'scaf' )
	out += ('' if isShort else 'fold' )
	out += ('_' if isUnSc else '' )
	if isEmpty:
		out = ''
	nname = name.lower().replace('scaffold', '').replace('scaf','').replace('_','')
	try:
		iname = int( nname )
		if isZero:
			out += '{number:0{width}d}'.format(width=(isZero+1), number=iname)
		else:
			out += str(iname)
	except ValueError:
		out += nname
		print( 'scaffold error number', nname )
	return out

def decodeContigOptions( options ):
	opts = list(options.lower())
	isCap = 'c' in opts
	isUnSc = 'u' in opts
	isZero = opts.count( 'z' )
	isEmpty = 'e' in opts
	return isCap, isUnSc, isZero, isEmpty

def formatContig( name, isCap, isUnSc, isZero, isEmpty ):
	out = ('Contig' if isCap else 'contig' )
	out += ('_' if isUnSc else '' )
	if isEmpty:
		out = ''
	nname = name.lower().replace('contig', '').replace('_','')
	try:
		iname = int( nname )
		if isZero:
			out += '{number:0{width}d}'.format(width=(isZero+1), number=iname)
		else:
			out += str(iname)
	except ValueError:
		out += nname
		print( 'contig error number', nname )
	return out

def decodeCLMOptions( options ):
	opts = list(options)
	mtType = None
	chType = None
	lmType = None
	for op in opts:
		# mitochondria
		if op.lower() in ['m','t','i']:
			if mtType != None:
				print( 'ERROR: too many specifications for mitochondria formatting' )
				exit()
			elif op == 'M':
				mtType = 'ChrM'
			elif op == 'm':
				mtType = 'chrM'
			elif op == 'T':
				mtType = 'ChrMt'
			elif op == 't':
				mtType = 'chrMt'
			elif op == 'I':
				mtType = 'Mitochondria'
			elif op == 'i':
				mtType = 'mitochondria'
		# chrloroplast
		elif op.lower() in ['c','h']:
			if chType != None:
				print( 'ERROR: too many specifications for chloroplast formatting' )
				exit()
			elif op == 'C':
				chType = 'ChrC'
			elif op == 'c':
				chType = 'chrC'
			elif op == 'H':
				chType = 'Chloroplast'
			elif op == 'h':
				chType = 'chloroplast'
		#lambda
		elif op.lower() in ['l','a']:
			if lmType != None:
				print( 'ERROR: too many specifications for chloroplast formatting' )
				exit()
			elif op == 'L':
				lmType = 'ChrL'
			elif op == 'l':
				lmType = 'chrL'
			elif op == 'A':
				lmType = 'Lambda'
			elif op == 'a':
				lmType = 'lambda'
	# end for
	return mtType, chType, lmType

def formatCLM( name, mtType, chType, lmType ):
	n=name.lower()
	if n in ['chrc', 'chloroplast','c']:
		return chType
	elif n in ['chrm','mitochrondria','chrmt','m']:
		return mtType
	elif n in ['chrl','lambda','l']:
		return lmType

def decodeOtherOptions( options ):
	opts = list(options.lower())
	
	isCap = 'c' in opts
	isLower = 'l' in opts
	if 'c' in opts and 'l' in opts:
		print( 'WARNING: do not specify capitalizing and lowercase for other..using neither' )
		isCap = False
		isLower = False
	addChr = 'h' in opts
	return isCap, isLower, addChr

def formatOther( name, isCap, isLower, addChr ):
	
	if addChr and isLower:
		out = 'chr'
	else:
		out = ('Chr' if addChr else '')
	if isCap:
		out += name[0].upper() + name[1:]
	elif isLower and not addChr:
		out += name[0].lower() + name[1:]
	else:
		out += name
	return out
	
