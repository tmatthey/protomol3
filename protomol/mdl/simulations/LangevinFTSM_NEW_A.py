# FTSM simulation of alanine dipeptide

from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from FTSM import *

# PHYSICAL SYSTEM
from Numeric import *

#import mpi
import time
import sys

#if (mpi.rank == 0):
#   start = time.time()
# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "examples/diAlanine/alanC5axial.pdb")
io.readPSF(phys, "examples/diAlanine/blockdialanine.psf")
io.readPAR(phys, "examples/diAlanine/par_all27_prot_lipid.par")
phys.bc = "Vacuum"
phys.temperature = 272
phys.exclude = "scaled1-4"
phys.cellsize = 6.5

forces = Forces()
ff = forces.makeForceField(phys, "charmm")
ff.constrainDihedral(11, 18)
    
# NOW RUN FTSM
iter = 0
numpoints = 18 #40
prop = Propagator()
ftsm = FTSM(phys, forces, prop, io, ff, 0.3)  # 40000 * 0.2 fs = 2 ps

# DETERMINE STARTING AND ENDING POINTS OF THE STRING
#workpt = mpi.rank

# INITIALIZE FTSM ON ROOT
print "INIT"
#if (mpi.rank == 0):
S = ftsm.initialize("examples/diAlanine/alanC5axial.pdb",
                    "examples/diAlanine/alanC7axial.pdb", numpoints)
print "\nI"+str(iter+1)+": ",
for jj in range(0, S.__len__()):
         print str(S[jj][0])+","+str(S[jj][1])+",",
sys.stdout.flush()

io.readPDBPos(phys, "examples/diAlanine/alanC5axial.pdb")
print "RUN"
# RUN 45 iterations of FTSM
while (iter < 1000):

    # WALK AND MAKE INITIAL STATES
    #if (mpi.rank == 0 and iter == 0):
       # FUTURE: CHECK UNITS OF C2
       # MAKE SURE RESTRAINT IS NOT STIFF
       # BE SURE IT ADEQUETLY COMPENSATES FOR GRADIENT OF V.
       # FUTURE: REMOVE TUBE CONSTRAINT
       # RIGHT NOW: SET U VERY HIGH, POTENTIAL WILL ALWAYS BE 0
    #   print "WALKING"
    #   ftsm.walk(S, c2=10, u=100, v=1.8)
    #   print "DONE"
    #mpi.barrier()
    
    # ON EACH PROCESSOR, UPDATE X
    # FUTURE: REMOVE TUBE CONSTRAINT
    # RIGHT NOW: SET U VERY HIGH, POTENTIAL WILL ALWAYS BE 0
    # FUTURE: REMOVE SLOPES (NOW SET TO 1 AND NOT USED)
    #if (iter < 50):
    #   C2 = 10
    #else:
    #   C2 = 500
    ftsm.run(S[15], 0, 1, 1, c1=0, c2=1000,u=100, v=1.8)
    
    # NOW ON EACH PROCESSOR, UPDATE Z
    #newpair = ftsm.changeZ(S[workpt], c2=100, alpha=200)
    
    # CREATE LIST ON THIS NODE, PASS BACK TO ROOT
    #currS = [newpair]
    #tempS = mpi.allgather(currS)
    
    # REPARAM ON ROOT
    # FUTURE DO NOT RETURN SLOPES
    #if (mpi.rank == 0):
    #   tS = ftsm.arclen(tempS, numpoints+2)
    #else:
    #   tS = []
    
    #S = mpi.bcast(tS)

    # OUTPUT S
    #if (mpi.rank == 0):
    #   print S.__len__()
    #   print "\nU"+str(iter+1)+": ",
    #   print S
    #   sys.stdout.flush()
     
    print iter,
    print ": ",
    print ftsm.ff.getPhiPsi() 
    iter += 1
    #mpi.barrier()

#mpi.synchronizedWriteln(str(mpi.rank)+": "+str(ftsm.ff.getPhiPsi()))
# TIME OF SIMULATION
#if (mpi.rank == 0):
#   stop = time.time()
#   print "TOTAL TIME: ",
#   print stop-start
