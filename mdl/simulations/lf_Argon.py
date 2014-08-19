# USING THE NEW STRUCTURE
from MDL import *

#from propagators.examplesP.leapfrog import *


# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/argon_280/untitled.pdb")
io.readPSF(phys, "data/argon_280/argon.psf")
io.readPAR(phys, "data/argon_280/argon.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.exclude = "scaled1-4"
phys.temperature = 106


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.nonbondedForces("l")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C1',
                             'cutoff':8.0}

# OUTPUT
#io.plots = {'totalenergy':2}
io.files = {'energies':('argon.energies', 1),
            'dcdtrajpos':('argon.pos.dcd', 2)}
io.screen = 1

# EXECUTE
prop = Propagator(phys, forces, io)
t0 = time.time()
#gamma = prop.propagate(scheme="leapfrog", steps=4000, dt=1.0, forcefield=ff, gpu=True)
gamma = prop.propagate(scheme="leapfrog", steps=4, dt=1.0, forcefield=ff, gpu=True)
t1 = time.time()
print "Execution Time in Seconds: ", t1 - t0
