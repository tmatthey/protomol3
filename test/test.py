import sys,os,glob,shlex,subprocess,comparator, argparse

parser = argparse.ArgumentParser(description='ProtoMol Test Suite')
parser.add_argument( '--verbose', '-v', action='store_true', default = False, help = 'Verbose output' )

group = parser.add_mutually_exclusive_group()
group.add_argument( '--single', '-s', help = 'Single test to run. Must be within the tests directory.' )
group.add_argument( '--regex', '-r', help = 'Regular expression of tests to run, Requires quotation marks around argument' )

args = parser.parse_args()

files = []

if args.single == None and args.regex == None:
	files=glob.glob("tests/*.conf")
else:
	if args.single != None:
		if args.single.find( ".conf" ) != -1:
			files.append( args.single )
		else:
			print("Invalid test configuartion")
			sys.exit(1)
	
	if args.regex != None:
		files=glob.glob("tests/"+args.regex+".conf")
		
files.sort()

print( files )

pwd = os.getcwd()

tests = 0
testspassed = 0
testsfailed = 0
failedtests = []

def printv( message ):
	if args.verbose == True:
		print( message )

for file in files:
	base = os.path.splitext( os.path.basename( file ) )[0]

	if os.path.exists( os.path.join( pwd, "ProtoMol" ) ) == True:
		print( "Executing Test: " + base )
		p = subprocess.Popen( shlex.split( os.path.join( pwd, "ProtoMol" ) + " " + file ), stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
	else:
		if( os.path.exists( os.path.join( pwd, "ProtoMol.exe" ) ) ) == True:
			print( "Executing Test: " + base )
			cmd = os.path.join( pwd, "ProtoMol.exe" ) + " " + file
			p = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
		else:
			print( "ProtoMol missing. Please put the ProtoMol executable in this directory" )
			sys.exit(1)
			
	expects = []
	outputs = glob.glob( "tests/output/" + base + ".*" )
	outputtemp = []
	for output in outputs:
		outbase = os.path.basename( output )
		if os.path.exists( "tests/expected/" + outbase ) == True:
			outputtemp.append( output )
			expects.append( "tests/expected/" + outbase )

	outputs = outputtemp

	for i in range( len( outputs ) ):
		tests += 1
		ftype = os.path.splitext( os.path.basename( outputs[i] ) )[1] 

		if ftype != ".dcd" and ftype != ".header" : 
			if ftype == ".vec":
				printv( "Testing: " + expects[i] + " " + outputs[i] )
			
				if comparator.compare( expects[i], outputs[i], "0.00001", 1.0, True, verbose = args.verbose ):
					printv( "Passed" )
					testspassed += 1
				else:
					printv( "Failed" )
					testsfailed += 1
					failedtests.append( "Comparison of " + expects[i] + " and " + outputs[i] )
			else:
				printv( "Testing: " + expects[i] + " " + outputs[i] )

				if comparator.compare( expects[i], outputs[i], "0.00001", verbose = args.verbose ):
					printv( "Passed" )
					testspassed += 1
				else:
					printv( "Failed" )
					testsfailed += 1
					failedtests.append( "Comparison of " + expects[i] + " and " + outputs[i] )

testsnotrun = tests - (testspassed + testsfailed )

print( "" )
print( "Tests: %d" % ( tests ) )
print( "Tests Not Run: %d" % ( testsnotrun ) )
print( "" )
print( "Tests Passed: %d" % ( testspassed ) )
print( "Tests Failed: %d\n" % ( testsfailed ) )

if len( failedtests ) > 0:
	print("")

	for i in range( len( failedtests ) ):
		print( "%s failed." % ( failedtests[i] ) )
