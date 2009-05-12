# Recursive Multiple Thermostat technique
from MDL import *

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/bpti_water_1101/bpti.pdb")
io.readPSF(phys, "data/bpti_water_1101/bpti.psf")
io.readPAR(phys, "data/bpti_water_1101/bpti.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.gbsa = True   # Generalized Born

# OUTPUT
io.screen = 2

# PROPAGATION
prop = Propagator(phys, forces, io)
prop.propagate(scheme="OpenMM", steps=200, dt=0.5, forcefield=ff)



