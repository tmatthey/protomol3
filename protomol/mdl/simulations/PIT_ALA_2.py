# A DRAFT OF A SIMULATION OF BLOCKED ALANINE DIPEPTIDE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from numpy import * # to generate random numbers
from pylab import semilogy,clf # to use matplotlib

def converge(xf,xc,no,i,n):
	#STRONG CONVERGENCE CRITERION TEST - TRAJECTORY DIFFERENCE
	nn = []
	for i in range(i,n+1):
		nn.append(linalg.norm(xf[i].positions-xc[i].positions))
	no.append(nn)

# DEFINE THE NUMBER OF COARSE POINTS
numpoints = 16
# DEFINE THE NUMBER OF ITERATIONS
numiter = 4

finestep = 0.025          # in fs
coarse2fine = 4          # ratio
coarsestep = finestep * coarse2fine #in fs
fixcoarse = 21
fixfine = 10
gammafine = 80
gammacoarse = gammafine
fdof = 0


# PHYSICAL SYSTEMS
io = IO()
normary = [] 
phif = []

xc0 = []
xf = []
dx = []
xc = []
forceary = []
propary = []
# 2 forces - coarse and fine
# 2 propagators - coarse and fine
forceary.append(Forces())
forceary.append(Forces())

phys = Physical()
io.readPDBPos(phys, "data/ALA/minC7eq.pdb")
io.readPSF(phys, "data/ALA/alan_mineq.psf")
io.readPAR(phys, "data/ALA/par_all27_prot_lipid.inp")
io.readEigenvectors(phys, "data/ALA/eigVmC7eq")
phys.bc = "Vacuum"
phys.temperature = 300
phys.exclude = "scaled1-4"
phys.seed = 1234

io = IO()	
for i in range(0, numpoints+1):
    xc0.append(phys.copy())
    dx.append(phys.copy())
    xc.append(phys.copy())
    xf.append(phys.copy())

dof = xc0[0].numAtoms()*3
temp = numpy.ndarray(dof)
tempv = numpy.ndarray(dof)
temp2 =numpy.ndarray(dof)
temp2v =numpy.ndarray(dof)

ff = forceary[0].makeForceField(xc0[0])
ff.bondedForces("badi")
ff.nonbondedForces("le")
ff.params = {'LennardJones':{'algorithm':'cutoff',
			     'switching':'C2',
			     'cutoff':12,
			     'switchon':9},
	     'CoulombDiElec':{
	
ff.setSwitching("LennardJones", "C2")
ff.setAlgorithm("LennardJones", "Cutoff")
ff.setCutoff("LennardJones", 12);
ff.setParameters("LennardJones", "switchon=9.0")
ff.setSwitching("CoulombDiElec", "C2")
ff.setAlgorithm("CoulombDiElec", "Cutoff")
ff.setCutoff("CoulombDiElec", 12);
ff.setParameters("CoulombDiElec", "switchon=9.0")
ff2 = forceary[0].makeForceField(xc0[0])
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

finef = forceary[1].makeForceField(dx[0])
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
finef2 = forceary[1].makeForceField(dx[0])
finef2.bondedForces("badi")
finef2.nonbondedForces("le")
finef2.setSwitching("LennardJones", "C2")
finef2.setAlgorithm("LennardJones", "Cutoff")
finef2.setCutoff("LennardJones", 5.5);
finef2.setParameters("LennardJones", "switchon=4.5")
finef2.setSwitching("CoulombDiElec", "C2")
finef2.setAlgorithm("CoulombDiElec", "Cutoff")
finef2.setCutoff("CoulombDiElec", 5.5);
finef2.setParameters("CoulombDiElec", "switchon=4.5")

#io.printScreen(1)	

#one step of PIT
#DEBUG
#for i in range(0,len(xc0)):
#	print "norm xc0[%d] %f" % (i,linalg.norm(xc0[i].positions))
#	print "norm xc[%d] %f" % (i,linalg.norm(xc[i].positions))
#	print "norm dx[%d] %f" % (i,linalg.norm(dx[i].positions))
#END DEBUG




#equilibrate
xf[i].seed=int(round(numpy.random.uniform(0,100000)))
forceary[1].forcevec.clear()
propary.append(Propagator(xc0[i], forceary[0]))
propary.append(Propagator(xf[i], forceary[1]))
propary[1].propagate("NormModeInt", xf[i], forceary[1], io, coarse2fine*100, 1, finef, ('fixmodes', fixfine), ('gamma', gammafine), ('fdof', fdof), finestep, finef2, ('avModeMass', 3.0),('seed',1017),('dtratio',coarse2fine))
for ii in range(0,dof):
	xc0[0].positions[ii]=xf[0].positions[ii]
	xc0[0].velocities[ii]=xf[0].velocities[ii]
	xc[0].positions[ii]=xf[0].positions[ii]
	xc[0].velocities[ii]=xf[0].velocities[ii]

#create numpoints seeds for coarse and fine propagators
seeds=numpy.random.uniform(0,numpoints*2,numpoints)
for i in range(0,len(seeds)):
	seeds[i]=int(round(seeds[i]))
#initial step: line 1 in pseudocode
for i in range(0,numpoints):
#	print "past line 1: step %d" % (i)
#coarse steps saved into xc0
#line 2 in pseudocode
	for ii in range(0,dof):
		temp[ii]=xc0[i].positions[ii]
		tempv[ii]=xc0[i].velocities[ii]
	forceary[0].forcevec.clear()
	propary[0].propagate("NormModeInt", xc0[i], forceary[0], io, 1, 1, ff, ('fixmodes', fixcoarse), ('gamma', gammacoarse), ('fdof', fdof), coarsestep, ff2, ('avModeMass',3.0), ('seed',seeds[i]),('dtratio',coarse2fine))
	for ii in range(0, dof):
		xc0[i+1].positions[ii] = xc0[i].positions[ii]
		xc0[i+1].velocities[ii] = xc0[i].velocities[ii]
		xc0[i].positions[ii] = temp[ii]
		xc0[i].velocities[ii] = tempv[ii]
		temp2[ii]=xf[i].positions[ii]
		temp2v[ii]=xf[i].velocities[ii]
		xf[i].positions[ii] = temp[ii]
		xf[i].velocities[ii] = tempv[ii]
#line 5 in pseudocode
	forceary[1].forcevec.clear()
        propary[1].propagate("NormModeInt", xf[i], forceary[1], io, coarse2fine, 1, finef, ('fixmodes', fixfine), ('gamma', gammafine), ('fdof', fdof), finestep, finef2, ('avModeMass', 3.0),('seed',seeds[i]),('dtratio',1))

#	phif.append(gamma[0].copy())
   	for ii in range(0, dof):
		xf[i+1].positions[ii] = xf[i].positions[ii]
		xf[i+1].velocities[ii] = xf[i].velocities[ii]
        	dx[i+1].positions[ii] = xf[i+1].positions[ii]-xc0[i+1].positions[ii]
		dx[i+1].velocities[ii] = xf[i+1].velocities[ii]-xc0[i+1].velocities[ii]
		xf[i].positions[ii] = temp2[ii]
		xf[i].velocities[ii] = temp2v[ii]


#	print "past line 5: step %d" % (i)
#print "finished initialization"
#DEBUG
#for i in range(0,len(xc0)):
#	print "norm xc0[%d] %f" % (i,linalg.norm(xc0[i].positions))
#	print "norm xc[%d] %f" % (i,linalg.norm(xc[i].positions))
#	print "norm dx[%d] %f" % (i,linalg.norm(dx[i].positions))
#END DEBUG
converge(xf,xc0,normary,1,numpoints)
#main loop
for j in range(numiter):
#create numpoints seeds for coarse and fine propagators
#	seeds=numpy.random.uniform(0,numpoints*2,numpoints)
#	for i in range(0,len(seeds)):
#		seeds[i]=int(round(seeds[i]))
	#line 9 in pseudocode
	for i in range(j,numpoints):
#		print "past line 9: step %d" % (i)
		#line 10 in pseudocode
		for ii in range(0,dof):
			temp[ii]=xc[i].positions[ii]
			tempv[ii]=xc[i].velocities[ii]
		forceary[0].forcevec.clear()
		propary[0].propagate("NormModeInt", xc[i], forceary[0], io, 1, 1, ff, ('fixmodes', fixcoarse), ('gamma', gammacoarse), ('fdof', fdof), coarsestep, ff2, ('avModeMass',3.0), ('seed',seeds[i]),('dtratio',coarse2fine))
		for ii in range(0, dof):
			xc0[i+1].positions[ii] = xc[i].positions[ii]
			xc0[i+1].velocities[ii] = xc[i].velocities[ii]
			xc[i+1].positions[ii] = xc0[i+1].positions[ii]+dx[i+1].positions[ii]
			xc[i+1].velocities[ii] = xc0[i+1].velocities[ii]+dx[i+1].velocities[ii]
			xc[i].positions[ii] = temp[ii]
			xc[i].velocities[ii] = tempv[ii]
			temp2[ii]=xf[i].positions[ii]
			temp2v[ii]=xf[i].velocities[ii]
			xf[i].positions[ii] = temp[ii]
			xf[i].velocities[ii] = tempv[ii]
		forceary[1].forcevec.clear()
		propary[1].propagate("NormModeInt", xf[i], forceary[1], io, coarse2fine, 1, ff, ('fixmodes', fixfine), ('gamma', gammafine), ('fdof', fdof), finestep, ff2, ('avModeMass', 3.0),('seed',seeds[i]),('dtratio',1))
   		for ii in range(0, dof):
			xf[i+1].positions[ii]=xf[i].positions[ii]
			xf[i+1].velocities[ii]=xf[i].velocities[ii]
			dx[i+1].positions[ii]=xf[i+1].positions[ii]-xc0[i+1].positions[ii]
			dx[i+1].velocities[ii]=xf[i+1].velocities[ii]-xc0[i+1].velocities[ii]
			xf[i].positions[ii]=temp2[ii]
			xf[i].velocities[ii]=temp2v[ii]
#		print "past line 14: step %d" % (i)
	converge(xf,xc,normary,1,numpoints)		

		
#analysis of error
pylab.clf()
x=range(len(normary[0]))
n=normary
pylab.semilogy(x,n[0],x,n[1],x,n[2])
pylab.xlabel('trajectory slice')
pylab.ylabel('norm of position error')
pylab.legend(('0th it','1st it','2nd it'),loc='lower right')
pylab.title('Coarse & fine NML DT=2dt GammaC=GammaF/2')

