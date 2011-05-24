import comparator, sys, argparse

parser = argparse.ArgumentParser(description='File comparison tool for floating point data')
parser.add_argument( 'original', help = 'Original file to compare against' )
parser.add_argument( 'new', help = 'New file to comare against' )
parser.add_argument( 'epsilon', type = float, help = 'Epsilon to compare files with' )
parser.add_argument( '--scale', type = float, default = float(1), help = 'Scaling factor for New file values' )
parser.add_argument( '--ignore-sign', action='store_true', default = False, help = 'Ignore sign differences between files' )
parser.add_argument( '--verbose', '-v', action='store_true', default = False, help = 'Verbose output' )

args = parser.parse_args()

result = comparator.compare( args.original, args.new, args.epsilon, args.scale, args.ignore_sign, args.verbose )

if result == False:
	print( "Files Differ" )
