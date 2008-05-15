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
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

forces = Forces()
ff = forces.makeForceField(phys)
ff2 = forces.makeForceField(phys, "charmm")

io.plots = {'totalenergy':4}
io.screen = 1
io.files = {'momentum':('somefile', 4),
            'energies':('enefile', 2)}

prop = Propagator(phys, forces, io)
prop.propagate(scheme="ShadowHMC", steps=100, cyclelength=20, dt=1, forcefield=[ff, ff2])

    
