# A DRAFT OF A SIMULATION OF 4-ATOM BUTANE
# USING THE NEW STRUCTURE
from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

phys = Physical()
io = IO()
io.readPDBPos(phys, "data/ALA/minC7eq.pdb")
io.readPSF(phys, "data/ALA/alan_mineq.psf")
io.readPAR(phys, "data/ALA/par_all27_prot_lipid.inp")
io.readEigenvectors(phys, "data/ALA/eigVmC7eq")
phys.bc = "Vacuum"
phys.temperature = 300
phys.cellsize = 5.0
phys.exclude = "scaled1-4"
phys.seed = 1234

forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("lc")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'cutoff':12,
                             'switchon':9}
ff.params['Coulomb'] = {'algorithm':'SCPISM',
                        'switching':'Cutoff',
                        'cutoff':12}
#ff.params['CoulombDiElec'] = ff.params['LennardJones'].copy()  # SAME AS LJ

ff2 = forces.makeForceField(phys)
ff2.bondedForces("badi")
ff2.nonbondedForces("lc")

ff2.params['LennardJones'] = {'algorithm':'Cutoff',
                              'switching':'C2',
                              'cutoff':5.5,
                              'switchon':4.5}

ff.params['Coulomb'] = {'algorithm':'SCPISM',
                        'switching':'Cutoff',
                        'cutoff':5.5,
                        'bornswitch':1}


#ff2.params['CoulombDiElec'] = ff2.params['LennardJones'].copy()  # SAME AS LJ

#io.plots = {'pressure':4,
#            'kineticenergy':1,
#            'temperature':1}
io.screen = 1

prop = Propagator(phys, forces, io)
prop.propagate(scheme=["NormModeInt", "NormModeMin"], 
	       steps=50, cyclelength=1, dt=4.0, forcefield=[ff, ff2], 
	       params={'NormModeInt':{'fixmodes':40,'gamma':91,'fdof':0},
		       'NormModeMin':{'avModeMass':30}}
              )


