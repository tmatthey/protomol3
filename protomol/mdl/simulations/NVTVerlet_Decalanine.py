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
io.readPDBPos(phys, "data/decalanine_66/alanine.pdb")
io.readPSF(phys, "data/decalanine_66/alanine.psf")
io.readPAR(phys, "data/decalanine_66/alanine.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("c")
ff.params['Coulomb'] = {'algorithm':'Cutoff',
                        'switching':'C1',
                        'cutoff':8.0}

# OUTPUT
io.plots = {'totalenergy':4}
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="NVTVerlet", steps=20, dt=0.5, forcefield=ff)
    
