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
#io.doMPL = 1
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Vacuum"
phys.cellsize = 6.5
phys.temperature = 300

forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badih")
ff.nonbondedForces("lc")

#ff2 = forces.makeForceField(phys)
#ff2.bondedForces("badih")
#ff2.nonbondedForces("lc")

ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
                                    'switching':['C2', 'C1'],
                                    'cutoff':8.0}

#hd = HDForce(phys,forces,3.0,1,1)
#ff.addPythonForce(hd)
#ff2.addPythonForce(hd)
ff.params['HarmonicDihedral'] = {'kbias':1,
                                 'dihedralnum':0,
                                 'angle':3.0
                                 }

#ff2.params['HarmonicDihedral'] = ff.params['HarmonicDihedral'].copy()

#io.plots = {'totalenergy':2}
io.files = {'xyztrajforce':('force.harm', 1)}
io.screen = 1

prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="Leapfrog", steps=40, dt=0.1, forcefield=ff)
