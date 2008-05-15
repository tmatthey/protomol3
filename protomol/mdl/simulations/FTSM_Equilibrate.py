from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

import FTSM
import numpy
io = IO()

# DEFINE THE NUMBER OF POINTS ON THE STRING 


# PHYSICAL SYSTEM
x = Physical()

io.readPDBPos(x, "data/alanDipeptideSol/solvate.pdb")
io.readPSF(x, "data/alanDipeptideSol/solvate.psf")
io.readPAR(x, "data/alanDipeptideSol/par_all27_prot_lipid.inp")
x.bc = "Periodic"
x.temperature = 300
x.exclude = "scaled1-4"
x.seed = 1234
dt = 1
T = 300


forces = Forces()
ff = forces.makeForceField(x)
ff.bondedForces("badi")
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



prop = Propagator(x, forces, io) 

print "EQUILIBRATING...50 ps"
prop.propagate("equilibrate", steps=50000, dt=dt, forcefield=ff, params={'T0':300, 'startatom':22})

# 500 ps
print "EQUILIBRATING...500 ps"
prop.propagate("equilibrate", steps=500000, dt=dt, forcefield=ff, params={'T0':300})


print "WRITING PDB"

io.writePDBPos(x, "data/alanDipeptideSol/solvate_eq.pdb")

print "DONE"
