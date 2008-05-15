# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/argon_280/untitled.pdb")
io.readPSF(phys, "data/argon_280/argon.psf")
io.readPAR(phys, "data/argon_280/argon.par")
#io.doMPL = True
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300
#phys.buildAll()

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.nonbondedForces("l")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                                    'switching':'C1',
                                    'cutoff':8.0}

# OUTPUT
io.plots = {'totalenergy':4}
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate(scheme="SelfconsistentLeapfrog", steps=200, dt=0.5, forcefield=ff)

