# A DRAFT OF A SIMULATION OF BLOCKED ALANINE DIPEPTIDE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from numpy import * # to generate random numbers
#from pylab import * # to use matplotlib

io = IO()
PHI_DIHEDRAL = 11
PSI_DIHEDRAL = 18

# DEFINE THE NUMBER OF COARSE POINTS
numpoints = 3
# DEFINE THE NUMBER OF ITERATIONS
numiter = 5

# PHYSICAL SYSTEMS
normarray = []
weakarray = []
physarray = []
finearray = []
resultarray = []
forcearray = []
proparray = []
# 2 forces - coarse and fine
# 2 propagators - coarse and fine
forcearray.append(Forces())
forcearray.append(Forces())
proparray.append(Propagator())
proparray.append(Propagator())
io = IO()	
for i in range(0, numpoints+1):
    physarray.append(Physical())
    finearray.append(Physical())
    resultarray.append(Physical())
    io.readPDBPos(physarray[i], "data/ALA/minC7eq.pdb")
    io.readPSF(physarray[i], "data/ALA/alan_mineq.psf")
    io.readPAR(physarray[i], "data/ALA/par_all27_prot_lipid.inp")
    io.readEigenvectors(physarray[i], "data/ALA/eigVmC7eq")
    physarray[i].bc = "Vacuum"
    physarray[i].temperature = 300
    physarray[i].exclude = "scaled1-4"
    physarray[i].seed = 1234+i
    io.readPDBPos(finearray[i], "data/ALA/minC7eq.pdb")
    io.readPSF(finearray[i], "data/ALA/alan_mineq.psf")
    io.readPAR(finearray[i], "data/ALA/par_all27_prot_lipid.inp")
    io.readEigenvectors(finearray[i], "data/ALA/eigVmC7eq")
    finearray[i].bc = "Vacuum"
    finearray[i].temperature = 300
    finearray[i].exclude = "scaled1-4"
    finearray[i].seed = 1234+i
    io.readPDBPos(resultarray[i], "data/ALA/minC7eq.pdb")
    io.readPSF(resultarray[i], "data/ALA/alan_mineq.psf")
    io.readPAR(resultarray[i], "data/ALA/par_all27_prot_lipid.inp")
    io.readEigenvectors(resultarray[i], "data/ALA/eigVmC7eq")
    resultarray[i].bc = "Vacuum"
    resultarray[i].temperature = 300
    resultarray[i].exclude = "scaled1-4"
    resultarray[i].seed = 1234+i

ff = forcearray[0].makeForceField(physarray[0])
ff.bondedForces("badi")
ff.nonbondedForces("le")
ff.setSwitching("LennardJones", "C2")
ff.setAlgorithm("LennardJones", "Cutoff")
ff.setCutoff("LennardJones", 12);
ff.setParameters("LennardJones", "switchon=9.0")
ff.setSwitching("CoulombDiElec", "C2")
ff.setAlgorithm("CoulombDiElec", "Cutoff")
ff.setCutoff("CoulombDiElec", 12);
ff.setParameters("CoulombDiElec", "switchon=9.0")
ff2 = forcearray[0].makeForceField(physarray[0])
ff2.bondedForces("badi")
ff2.nonbondedForces("le")
ff2.setSwitching("LennardJones", "C2")
ff2.setAlgorithm("LennardJones", "Cutoff")
ff2.setCutoff("LennardJones", 5.5);
ff2.setParameters("LennardJones", "switchon=4.5")
ff2.setSwitching("CoulombDiElec", "C2")
ff2.setAlgorithm("CoulombDiElec", "Cutoff")
ff2.setCutoff("CoulombDiElec", 5.5);
ff2.setParameters("CoulombDiElec", "switchon=4.5")

finef = forcearray[1].makeForceField(physarray[0])
finef.bondedForces("badi")
finef.nonbondedForces("le")
finef.setSwitching("LennardJones", "C2")
finef.setAlgorithm("LennardJones", "Cutoff")
finef.setCutoff("LennardJones", 12);
finef.setParameters("LennardJones", "switchon=9.0")
finef.setSwitching("CoulombDiElec", "C2")
finef.setAlgorithm("CoulombDiElec", "Cutoff")
finef.setCutoff("CoulombDiElec", 12);
finef.setParameters("CoulombDiElec", "switchon=9.0")

io.printScreen(1)	

#one step of PIT
coarsestep = 4	# in fs
finestep = 0.5  # in fs
coarse2fine = coarsestep / finestep
#create numpoints seeds for coarse and fine propagators
seeds=numpy.random.uniform(0,numpoints*2,numpoints)
for i in range(0,len(seeds)):
	seeds[i]=int(round(seeds[i]))
#initial step: line 1 in pseudocode
for i in range(0,numpoints):
#	print "past line 1: step %d" % (i)
#coarse steps saved into physarray
#line 2 in pseudocode
	gamma = proparray[0].propagate("NormModeInt", physarray[i], forcearray[0], io, 1, 1, ff, ('fixmodes', 40), ('gamma', 91), ('fdof', 0), coarsestep, ff2, ('avModeMass', 3.0),('seed',seeds[i]))
   	for ii in range(0, physarray[0].numAtoms()*3):
        	physarray[i+1].positions[ii] = gamma[0][ii]
#	print "past line 2: step %d" % (i)
#line 5 in pseudocode
	gamma = proparray[1].propagate("LangevinImpulse",physarray[i],forcearray[1],io,coarse2fine,finestep,finef,('seed',1000),('gamma',91))

   	for ii in range(0, physarray[0].numAtoms()*3):
        	finearray[i+1].positions[ii] = gamma[0][ii]-physarray[i+1].positions[ii]
#	print "past line 5: step %d" % (i)
#print "finished initialization"
#main loop
for j in range(5):
#create numpoints seeds for coarse and fine propagators
	seeds=numpy.random.uniform(0,numpoints*2,numpoints)
	for i in range(0,len(seeds)):
		seeds[i]=int(round(seeds[i]))
	#line 9 in pseudocode
	for i in range(0,numpoints):
#		print "past line 9: step %d" % (i)
		#line 10 in pseudocode
		gamma = proparray[0].propagate("NormModeInt", resultarray[i], forcearray[0], io, 1, 1, ff, ('fixmodes', 40), ('gamma', 91), ('fdof', 0), coarsestep, ff2, ('avModeMass', 3.0),('seed',seeds[i]))
   		for ii in range(0, physarray[0].numAtoms()*3):
        		physarray[i+1].positions[ii] = gamma[0][ii]
			resultarray[i+1].positions[ii] = gamma[0][ii]
#		print "past line 10: step %d" % (i)
		#line 11 in pseudocode
		resultarray[i+1].positions+=finearray[i+1].positions
#		print "past line 11: step %d" % (i)
		#line 14 in pseudocode
		gamma = proparray[1].propagate("LangevinImpulse",resultarray[i],forcearray[1],io,coarse2fine,finestep,finef,('seed',1000),('gamma',91))
   		for ii in range(0, physarray[0].numAtoms()*3):
        		finearray[i+1].positions[ii] = gamma[0][ii]-physarray[i+1].positions[ii]
#		print "past line 14: step %d" % (i)
	#CONVERGENCE CRITERION TEST
	#STRONG CONVERGENCE CRITERION - POSITIONS
	nn=linalg.norm(resultarray[i+1].positions-physarray[i+1].positions)
	print "trajectory difference: %f" % (nn)
	normarray.append(nn)
	#ANOTHER CONVERGENCE CRITERION - POTENTIAL ENERGY
	wn=forcearray[0].energies.potentialEnergy()-forcearray[1].energies.potentialEnergy()
	print "PE difference: %f" % (wn)
	weakarray.append(wn)


			
		
		
	
