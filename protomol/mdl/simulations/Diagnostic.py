# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

import FTSM
import sys

io = IO()
PHI_DIHEDRAL = 11
PSI_DIHEDRAL = 18

# DEFINE THE NUMBER OF POINTS ON THE STRING 
numpoints = 40


# PHYSICAL SYSTEMS
phys = Physical()
io.readPDBPos(phys, "data/diAlanine/alanC5axial.pdb")
io.readPSF(phys, "data/diAlanine/blockdialanine.psf")
io.readPAR(phys, "data/diAlanine/par_all27_prot_lipid.inp")
phys.bc = "Vacuum"
phys.temperature = 300
phys.exclude = "scaled1-4"
phys.seed = 1234
forces = Forces()


ff = forces.makeForceField(phys)
ff.bondedForces("badihh") # 2 harmonic dihedral forces
ff.nonbondedForces("lc")

prop = Propagator(phys, forces, io)
# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "data/diAlanine/alanC5axial.pdb")

#kappa = 1000.0
#gamma = 50000.0
kappa = 40.0
gamma = 2000.0

#io.initializePlot('string')
#io.pause=1
stringgraph=io.newGraph('Phi', 'Psi')

# PRINT INITIAL STRING I0
io.plotVector(prop,stringgraph,[[phys.phi(PHI_DIHEDRAL), phys.phi(PSI_DIHEDRAL)]], rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")

dt = 1.0

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI_DIHEDRAL-1, PSI_DIHEDRAL-1],
                                 'angle':[phys.phi(PHI_DIHEDRAL), phys.phi(PSI_DIHEDRAL)]}

for iter in range(0, 200000):
    prop.propagate(scheme="LangevinImpulse", steps=1, dt=dt, forcefield=ff)
    io.plotVector(prop,stringgraph,[[phys.phi(PHI_DIHEDRAL), phys.phi(PSI_DIHEDRAL)]], rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
    
