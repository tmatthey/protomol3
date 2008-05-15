# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

import FTSM
import math

io = IO()
PHI=11
PSI=18
kappa = 40

# PHYSICAL SYSTEM
x = Physical()

io.readPDBPos(x, "data/alanDipeptideSol/solvate_eq.pdb")
io.readPSF(x, "data/alanDipeptideSol/solvate.psf")
io.readPAR(x, "data/alanDipeptideSol/par_all27_prot_lipid.inp")
x.bc = "Periodic"
x.temperature = 300
x.exclude = "scaled1-4"
x.seed = 1234

force = Forces()
prop = Propagator(x, force, io)

ff = force.makeForceField(x)
ff.bondedForces("badihh")
ff.nonbondedForces("lc")

ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'Cutoff',
                             'cutoff':9}

ff.params['Coulomb'] = {'algorithm':'PME',
                        'switching':'Cutoff',
                        'interpolation':'BSpline',
                        'gridsize':32,
                        'order':4,
                        'cutoff':9}

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa, kappa],
                                 'dihedralnum':[PHI-1, PSI-1, 3],
                                 'angle':[-1.5, -1.5, 0.4]}

stringgraph=io.newGraph('Phi', 'Psi')

for ii in range(0, 500):
    prop.propagate(scheme="velocityscale", steps=1, dt=1, forcefield=ff, params={'T0':300})
    phi = x.angle(PHI)
    psi = x.angle(PSI)
    io.plotVector(prop,stringgraph,[[phi,psi]], rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])


