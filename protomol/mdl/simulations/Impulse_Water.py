# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *


phys = Physical()
io = IO()
io.readPDBPos(phys, "data/water_216/water216.new.pos.pdb")
io.readPDBVel(phys, "data/water_216/water216.new.vel.pdb")
io.readPSF(phys, "data/water_216/water216.psf")
io.readPAR(phys, "data/water_216/water216.par")
phys.bc = "Periodic"
phys.temperature = 300
phys.cellsize = 6.5

forces = Forces()
ff = forces.makeForceField(phys)
ff.nonbondedForces("lc")
ff2 = forces.makeForceField(phys)
ff2.bondedForces("badi")

ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
                                    'switching':['C2', 'C1'],
                                    'cutoff':8.0}

io.plots = {'pressure':5}
io.screen = 2


prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme=["impulse", "takahashi"], steps=100, cyclelength=5, dt=0.1, forcefield=[ff, ff2])
