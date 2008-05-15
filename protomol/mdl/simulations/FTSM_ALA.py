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
numpoints = 20


# PHYSICAL SYSTEMS
physarray = []
forcearray = []
proparray = []
for i in range(0, numpoints):
    physarray.append([Physical(), Physical()])
    forcearray.append([Forces(), Forces()])
    for j in range(0, 2):
        if (i > numpoints / 2):
            io.readPDBPos(physarray[i][j], "data/diAlanine/alanC5axial.pdb")
        else:
            io.readPDBPos(physarray[i][j], "data/diAlanine/alanC7axial.pdb")
        io.readPSF(physarray[i][j], "data/diAlanine/blockdialanine.psf")
        io.readPAR(physarray[i][j], "data/diAlanine/par_all27_prot_lipid.inp")
        physarray[i][j].bc = "Vacuum"
        physarray[i][j].temperature = 300
        physarray[i][j].exclude = "scaled1-4"
        physarray[i][j].seed = 1234
    proparray.append([Propagator(physarray[i][0],forcearray[i][0],io), Propagator(physarray[i][1],forcearray[i][1],io)])


ff = forcearray[0][0].makeForceField(physarray[0][0])
ff.bondedForces("badihh") # 2 harmonic dihedral forces
ff.nonbondedForces("lc")

###################################################################
# STEP 1: OBTAIN THE INITIAL STRING
S = [] # INITIALIZE TO EMPTY
# PROPAGATE ONCE SO OUR CONSTRAINING FORCE CAN GIVE US THE INITIAL
# PHI, PSI
phiI = physarray[0][0].phi(PHI_DIHEDRAL)
psiI = physarray[0][0].phi(PSI_DIHEDRAL)
print phiI, " ", psiI

# READ THE SECOND PDB
# PROPAGATE ONCE TO GET THE FINAL PHI,PSI
#io.readPDBPos(physarray[0][0], "data/diAlanine/alanC7axial.pdb")
phiF = physarray[numpoints-1][0].phi(PHI_DIHEDRAL)
psiF = physarray[numpoints-1][0].phi(PSI_DIHEDRAL)
print phiF, " ", psiF

# SEPARATE INTO EQUIDISTANT PHI AND PSI
dphi = (phiF - phiI) / (numpoints-1)
dpsi = (psiF - psiI) / (numpoints-1)
S.append([phiI, psiI])
newphi = phiI
newpsi = psiI
for ii in range(1, numpoints-1):
    newphi += dphi
    newpsi += dpsi
    S.append([newphi, newpsi])
S.append([phiF, psiF])

# PRINT INITIAL STRING I0
print "\nI0: ",
print S

# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "data/diAlanine/alanC5axial.pdb")

kappa = 40.0
gamma = 2000.0

#io.initializePlot('string')
#io.pause=1
stringgraph=io.newGraph('Phi', 'Psi')

# PRINT INITIAL STRING I0
print "\nI0: ",
print S
io.plotVector(proparray[0][0],stringgraph,S, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")

dt = 1.0

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI_DIHEDRAL-1, PSI_DIHEDRAL-1],
                                 'angle':[S[0][0], S[0][1]]}
        
for iter in range(0, 200000): # NUMBER OF FTSM ITERATIONS
    for workpt in range(0, numpoints): # LOOPING OVER POINTS
        # UPDATE FREE SPACE
        # USE FIRST SYSTEM TO GET M
        # USE SECOND SYSTEM TO OBTAIN PHI AND PSI DIFFERENCES
        # FROM TARGETS
        #if (iter >= 10000 and iter <= 100000):
        #    kappa += (25000.-1000.)/90000.
        M_phiphi = FTSM.M_phi(physarray[workpt][0], PHI_DIHEDRAL)
        M_psipsi = FTSM.M_phi(physarray[workpt][0], PSI_DIHEDRAL)
        M_phipsi = FTSM.M_phipsi(physarray[workpt][0], PHI_DIHEDRAL, PSI_DIHEDRAL)
        #print "DIFFERENCE: ", S[workpt][0]-physarray[workpt][1].phi(PHI_DIHEDRAL), " " , S[workpt][1]-physarray[workpt][1].phi(PSI_DIHEDRAL)
        S[workpt][0] -= (kappa/gamma)*dt*(M_phiphi*(S[workpt][0]-physarray[workpt][1].phi(PHI_DIHEDRAL)) + M_phipsi*(S[workpt][1] - physarray[workpt][1].phi(PSI_DIHEDRAL)))
        S[workpt][1] -= (kappa/gamma)*dt*(M_phipsi*(S[workpt][0]-physarray[workpt][1].phi(PHI_DIHEDRAL)) + M_psipsi*(S[workpt][1] - physarray[workpt][1].phi(PSI_DIHEDRAL)))

        # UPDATE CARTESIAN
        # Dr. Izaguirre: I have checked and this constraint
        # is correct.  The energy is harmonic, but the force (the gradient)
        # is not harmonic.  In fact it is exactly what is in the paper.

        proparray[workpt][0].propagate(scheme="LangevinImpulse", steps=1, dt=dt, forcefield=ff)

        proparray[workpt][1].propagate(scheme="LangevinImpulse", steps=1, dt=dt, forcefield=ff)

        # My own function which sets phi and psi for individual force objects
        # Saves performance since I only change 'angle', I don't want to define
        # all new force objects by changing params.
        if (workpt == numpoints-1):
           FTSM.setPhiPsi(ff, S[0][0], S[0][1])
        else:
           FTSM.setPhiPsi(ff, S[workpt+1][0], S[workpt+1][1])
        
    # REPARAMETERIZE S
    #S = FTSM.arclen(S, numpoints, 1000)
    #print "\nI"+str(iter+1)+": ",S
    print "\nI"+str(iter+1)+" BEFORE SMOOTH: ",S
    # REPARAMETERIZE S
    #S.reverse()
    #FTSM.switchPhiPsi(S)
    S = FTSM.arclenChris(FTSM.extractPhi(S), FTSM.extractPsi(S))
    #FTSM.switchPhiPsi(S)
    #S.reverse()
    print "\nI"+str(iter+1)+": ",S
    io.plotVector(proparray[0][0],stringgraph,S, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
