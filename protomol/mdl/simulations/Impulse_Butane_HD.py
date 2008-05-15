# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from forces.HDForce import * # PYTHON FORCE CLASS

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Vacuum"
phys.cellsize = 6.5
phys.temperature = 300

forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("lc")

#ff2 = forces.makeForceField(phys)
#ff2.bondedForces("badi")
#ff2.nonbondedForces("lc")

#ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
#                                    'switching':['C2', 'C1'],
#                                    'cutoff':8.0}

hd = HDForce(phys,forces,3.0,1,1)  # THIS IS A PYTHON-PROTOTYPED FORCE
ff.addPythonForce(hd)

io.screen = 1
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="Leapfrog", steps=40, dt=0.1, forcefield=ff)
