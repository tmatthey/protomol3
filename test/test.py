import sys,os,glob,shlex,subprocess,comparator

files=glob.glob("tests/*.conf")

pwd = os.getcwd()

tests = 0
testspassed = 0
testsfailed = 0

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
				print( "Testing: " + expects[i] + " " + outputs[i] )
				if comparator.compare( expects[i], outputs[i], "0.00001", 1.0, True ):
					print( "Passed" )
					testspassed += 1
				else:
					print( "Failed" )
					testsfailed += 1
			else:
				print( "Testing: " + expects[i] + " " + outputs[i] )
				if comparator.compare( expects[i], outputs[i], "0.00001" ):
					print( "Passed" )
					testspassed += 1
				else:
					print( "Failed" )
					testsfailed += 1

print( "" )
print( "Tests: %d" % ( tests ) )
print( "Tests Passed: %d" % ( testspassed ) )
print( "Tests Failed: %d" % ( testsfailed ) )

testsnotrun = tests - (testspassed + testsfailed )
print( "Tests Not Run: %d" % ( testsnotrun ) )