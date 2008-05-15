from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *
from numpy import *


# PHYSICAL SYSTEM
phys = Physical()
phys2 = Physical()
io = IO()
io.readPDBPos(phys, "examples/UA_butane/UA_butane.pdb")
io.readPSF(phys, "examples/UA_butane/UA_butane.psf")
io.readPAR(phys, "examples/UA_butane/UA_butane.par")
io.readPDBPos(phys2, "examples/UA_butane/UA_butane.pdb")
io.readPSF(phys2, "examples/UA_butane/UA_butane.psf")
io.readPAR(phys2, "examples/UA_butane/UA_butane.par")
phys.bc = "Vacuum"
phys.exclude = "scaled1-4"
phys.cellsize = 6.5
phys.temperature = 300
phys2.bc = "Vacuum"
phys2.exclude = "scaled1-4"
phys2.cellsize = 6.5
phys2.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")
ff.setCutoff("LennardJones", 10.0)
ff.setSwitching("LennardJones", "C1")
ff.setAlgorithm("LennardJones", "Cutoff")
ff.setParameters("LennardJones", "blocksize=32")

# OUTPUT
#io.plotTemperatures(4)
#io.printScreen(10)
io.writeTemperature("output.temperatures", 2)
io.writeEnergies("output.energies", 4)
io.writeMomentum("output.momentum", 5)
io.writeDiffusion("output.diffusion", 3)
io.writePDBFramePos(phys, "output.frame.pdb", 2)
io.writeXYZTrajectoryForce("output.force.xyz", 2)
io.writeXYZTrajectoryForce("output.pos.xyz", 2)
io.writeXYZTrajectoryForce("output.vel.xyz", 2)

# PROPAGATION
prop = Propagator()
dof = phys.numAtoms()*3
temp = numpy.ndarray(dof)
out1 = numpy.ndarray(dof)
out2 = numpy.ndarray(dof)
for ii in range(0,dof):
    temp[ii]=phys.positions[ii]

phys.seed=100
prop.propagate("BBK", phys, forces, io, 1, 0.5, ff,('seed',100))

#print phys.positions
#for ii in range(0,dof):
#    out1[ii]=phys.positions[ii]

io.writePAR(phys, "output.par")
io.writePDBPos(phys, "final.pdb.pos")
io.writePDBVel(phys, "final.pdb.vel")
io.writeXYZPos(phys, "final.xyz.pos")
io.writeXYZVel(phys, "final.xyz.vel")
io.writeXYZBinPos(phys, "final.xyz.bin.pos")
io.writeXYZBinVel(phys, "final.xyz.bin.vel")

for ii in range(0,dof):
    phys2.positions[ii]=temp[ii]

prop1 = Propagator()

phys2.seed=100
forces.forcevec.clear()
prop1.propagate("BBK", phys2, forces, io, 1, 0.5, ff,('seed',100))

#print phys2.positions
#for ii in range(0,dof):
#    out2[ii]=phys2.positions[ii]

