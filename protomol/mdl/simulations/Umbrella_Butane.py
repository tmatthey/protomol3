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

ff2 = forces.makeForceField(phys)
ff2.bondedForces("badih")
ff2.nonbondedForces("lc")

ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
                                    'switching':['Universal', 'Universal'],
                                    'cutoff':8.0}
ff.params['HarmonicDihedral'] = {'kbias':1,
                                 'dihedralnum':0,
                                 'angle':355
                                 }

ff2.params['HarmonicDihedral'] = ff.params['HarmonicDihedral'].copy()

io.plots = {'totalenergy':2}
io.screen = 2

prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="Umbrella", steps=40, cyclelength=5, dt=0.1, forcefield=[ff, ff2])
