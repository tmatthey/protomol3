# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *


# TIMING THIS SIMULATION
import time
start = time.time()

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/bpti_water_1101/bpti.pdb")
io.readPSF(phys, "data/bpti_water_1101/bpti.psf")
io.readPAR(phys, "data/bpti_water_1101/bpti.par")
phys.bc = "Periodic"
phys.cellsize = 5
phys.exclude = "scaled1-4"
phys.temperature = 0
phys.seed = 7536031

# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")
ff.params['LennardJones'] = {'algorithm':'SimpleFull'}
ff.params['Coulomb'] = {'algorithm':'MultiGrid',
                        'switching':'C3'}

#ff.fastelectro['MultiGrid'] = {'interpolation':'Hermite',
#                               'toplevelgrid':8,
#                               'levels':2,
#                               's':10,
#                               'order':6,
#                               'ratio':2,
#                               'h':3,
#                               'origin':[0,0,0]}
#ff.fastelectro['MultiGrid'] = {'gridsize':32,
#                               'cutoff':16}

# OUTPUT
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="Leapfrog", steps=2000, dt=0.5, forcefield=ff)

# PRINT ELAPSED TIME
stop = time.time()
print "ELAPSED: ",
print stop-start
