# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE

from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from forces.HDForce import * # PYTHON FORCE CLASS
from forces.ElectrostaticForce import *
import time
start = time.time()

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanDipeptideBlock/alanC7axial.pdb")
io.readPSF(phys, "data/alanDipeptideBlock/blockdialanine.psf")
io.readPAR(phys, "data/alanDipeptideBlock/par_all27_prot_lipid.inp")
phys.bc = "Vacuum"
phys.temperature = 300
phys.exclude = "scaled1-4"
phys.cellsize = 226.5


# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")

# IO
io.screen = 1
io.plots = {'kineticenergy':2}
# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate(scheme="LangevinImpulse", steps=1000, dt=1.0, forcefield=ff)
    
stop=time.time()
print "TOTAL TIME: ", stop-start
