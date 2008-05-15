# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

phys = Physical()
io = IO()
io.readPDBPos(phys, "data/WW/wwd.min.pdb")
io.readPSF(phys, "data/WW/wwd_mineq.psf")
io.readPAR(phys, "data/WW/par_all27_prot_lipid.inp")
io.readEigenvectors(phys, "data/WW/wwdaevect.bin.dat")
phys.bc = "Vacuum"
phys.temperature = 300
phys.cellsize = 4.0
phys.exclude = "scaled1-4"
phys.seed = 11

forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("le")

ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'Cn',
                             'cutoff':10,
                             'switchon':8,
                             'order':4}
ff.params['CoulombDiElec'] = {'algorithm':'Cutoff',
                              'switching':'Cn',
                              'cutoff':10,
                              'switchon':0,
                              'order':4}

ff2 = forces.makeForceField(phys)
ff2.bondedForces("badi")
ff2.nonbondedForces("le")

ff2.params['LennardJones'] = {'algorithm':'Cutoff',
                              'switching':'Cn',
                              'cutoff':5,
                              'switchon':3.5,
                              'order':4}
ff2.params['CoulombDiElec'] = {'algorithm':'Cutoff',
                               'switching':'Cn',
                               'cutoff':10,
                               'switchon':0,
                               'order':4}

io.screen = 2

prop = Propagator(phys,forces,io)
prop.propagate(scheme=["NormModeInt", "NormModeMin"], 
	       steps=500, cyclelength=1, dt=10.0, forcefield=[ff, ff2], 
	       params={'NormModeInt':{'fixmodes':1400}}
              )
