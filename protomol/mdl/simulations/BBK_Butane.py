from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

import time
start = time.time()

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
#io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
#io.readPSF(phys, "data/UA_butane/UA_butane.psf")
#io.readPAR(phys, "data/UA_butane/UA_butane.par")
print "IO"
io.readPDBPos(phys, "data/bpti_water_1101/bpti.pdb")
print "PDB"
io.readPSF(phys, "data/bpti_water_1101/bpti.psf")
print "PSF"
io.readPAR(phys, "data/bpti_water_1101/bpti.par")
print "PAR"
phys.bc = "Periodic"
phys.exclude = "scaled1-4"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff2 = forces.makeForceField(phys)
ff.bondedForces("badi")
ff2.nonbondedForces("lc")
#ff = forces.makeForceField(phys, "charmm")
#ff.params['LennardJones'] = {'algorithm':'Cutoff',
#                             'switching':'C1',
#                             'cutoff':8.0}

# OUTPUT
#io.plots = {'totalenergy':4}
io.screen = 1
#io.files = {'temperature':("output.temperatures", 2),
#            'energies':("output.energies", 4),
#            'momentum':("output.momentum", 5),
#            'diffusion':("output.diffusion", 3),
#            'pdbframepos':("output.frame.pdb",2),
#            'xyztrajforce':("output.force.xyz",2,0),
#            'xyztrajpos':("output.pos.xyz",2),
#            'xyztrajvel':("output.vel.xyz",2)}

# PROPAGATION
prop = Propagator(phys, forces, io)
#gamma = prop.propagate(["impulse", "leapfrog"], 500, 4, 0.5, [ff2, ff])
#gamma = prop.propagate(scheme="Impulse", steps=1000, cyclelength=4, dt=0.5, forcefield=[ff2, ff])
gamma = prop.propagate(scheme="BBK", steps=100, dt=1.0, forcefield=ff, params={'shake':'on'})
stop=time.time()
print "TOTAL TIME: ", stop-start
#io.writePAR(phys, "output.par")
#io.writePDBPos(phys, "final.pdb.pos")
#io.writePDBVel(phys, "final.pdb.vel")
#io.writeXYZPos(phys, "final.xyz.pos")
#io.writeXYZVel(phys, "final.xyz.vel")
#io.writeXYZBinPos(phys, "final.xyz.bin.pos")
#io.writeXYZBinVel(phys, "final.xyz.bin.vel")
