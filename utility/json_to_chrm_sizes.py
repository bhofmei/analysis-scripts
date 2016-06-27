import sys, os, json
from io import open

# Usage: python json_to_chrm_sizes.py <in_json> <out_name>

def processInputs( inFileStr, outFileStr ):
	print( 'Converting {:s} to {:s}'.format( os.path.basename( inFileStr ), outFileStr ) )
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	jsonAr = json.load( inFile )
	for chrm in jsonAr:
		outStr = u'{:s}\t{:d}\n'.format( chrm['name'], chrm['length'] )
		outFile.write( outStr )
	inFile.close()
	outFile.close()
	

def parseInputs( argv ):
	inFileStr = argv[0]
	outFileStr = argv[1]
	if outFileStr.endswith('.tsv' ) == False:
		outFileStr += '.tsv'
	processInputs( inFileStr, outFileStr )


if __name__ == "__main__":
	if len(sys.argv) < 3 :
		print ("Usage: python json_to_fai.py <in_json> <out_name>\nif out_name doesn't end in .tsv, it will be added")
	else:
		parseInputs( sys.argv[1:] )
