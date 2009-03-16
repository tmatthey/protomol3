# Recursive Multiple Thermostat technique
from MDL import *

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
#io.doMPL = 1
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
print "FF"
#ff = forces.makeForceField(phys)
ff = forces.makeForceField(phys, "charmm")
print "GG"

# OUTPUT
#io.plots = {'kineticenergy':4}
print "HH"
io.screen = 2

print "II"
# PROPAGATION
print "SS"
prop = Propagator(phys, forces, io)
print "PROPAGATE"
prop.propagate(scheme="RMT", steps=200, dt=0.5, forcefield=ff)


    
