import comparator, sys

if len( sys.argv ) < 4:
	print( "Invalid number of arguments" )
	sys.exit( 1 )

epsilon = float( sys.argv[3] )

if len( sys.argv) == 5:
	comparator.compare( sys.argv[1], sys.argv[2], epsilon, float( sys.argv[4] ) )
else:
	comparator.compare( sys.argv[1], sys.argv[2], epsilon )
