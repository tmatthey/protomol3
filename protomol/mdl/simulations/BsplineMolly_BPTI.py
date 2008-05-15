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
io.readPDBPos(phys, "data/bpti_water_1101/bpti.pdb")
io.readPSF(phys, "data/bpti_water_1101/bpti.psf")
io.readPAR(phys, "data/bpti_water_1101/bpti.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

forces = Forces()
ff = forces.makeForceField(phys)
ff2 = forces.makeForceField(phys) 
ff.nonbondedForces("lc")
ff2.bondedForces("badi")
ff2.mollyForces("ba")

ff.params['LennardJonesCoulomb'] = {'algorithm':'Cutoff',
                                    'switching':['C2', 'C1'],
                                    'cutoff':8.0}
 
io.screen = 2

prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme=["BsplineMolly","BBK"], steps=20, cyclelength=5, dt=0.1, forcefield=[ff, ff2])

    
