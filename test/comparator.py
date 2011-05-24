import sys,math

def absDiff( one, two ):
	return math.fabs( one - two )
	
def sabsDiff( one, two ):
	return math.fabs( one ) - math.fabs( two )

def compare( fExpected, fNew, epsilon, scalar=1.0, ignoreSign = False, verbose = False ):
	fone = open( fExpected, 'r' )
	ftwo = open( fNew, 'r' )

	epsilon = float( epsilon )

	diffs = 0
	
	i = 0
	for lineone in fone:
		linetwo = ftwo.readline()

		elementsone = lineone.split()
		elementstwo = linetwo.split()

		i = i + 1

		lenone = len( elementsone )
		lentwo = len( elementstwo )

		if lenone != lentwo:
			diffs = diffs + 1
			if verbose == True:
				print( "Line: %d differs in size." % ( i ) )
				print( "Should be %d elements but there are %d." % ( lenone, lentwo ) )
		else:
			for j in range( lenone ):
				try:
					feone = float( elementsone[j] ) * scalar
					fetwo = float( elementstwo[j] )
					fediff = 0
					
					if ignoreSign == True:
						fediff = sabsDiff( feone, fetwo )
					else:
						fediff = absDiff( feone, fetwo )

					if fediff > epsilon:
						diffs = diffs + 1
						if verbose == True:
							print( "Line %d, Element %d Differs" % ( i, j ) )
							print( "Expected: %f, Actual: %f, Difference: %f" % ( feone, fetwo, fediff ) )
				except ValueError:
					x = 1
	return diffs == 0


