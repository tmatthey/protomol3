# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanine_gromacs/alanine.pdb")
io.readGromacs(phys, "data/alanine_gromacs/alanine.top", "data/alanine_gromacs/ffamber96/")
phys.bc = "Vacuum"
phys.temperature = 106


# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")

# OUTPUT
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("Leapfrog", steps=1000, dt=0.5, forcefield=ff)
