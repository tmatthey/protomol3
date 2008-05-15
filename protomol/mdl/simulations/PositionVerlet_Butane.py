# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE

from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Periodic"
phys.cellsize = 0.1
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")

ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
                                    'switching':['C2', 'C1'],
                                    'cutoff':8.0}

# OUTPUT
io.plots = {'totalenergy':4}
io.screen = 2

# PROPAGATION
prop = Propagator(phys, forces, io)
prop.propagate(scheme="PositionVerlet", steps=20, dt=0.5, forcefield=ff)
    
