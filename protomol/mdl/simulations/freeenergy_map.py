from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *


import sys
import os

def setPhiPsi(phi,psi,forcefield):
    flag = False
    for i in range(0, forcefield.forcetypes.__len__()):
        if (forcefield.forcetypes[i] == 'h'):
            if (not flag):
               forcefield.forcearray[i].setReferenceDihedral(phi)
               flag = True
            else:
               forcefield.forcearray[i].setReferenceDihedral(psi)

   

numgridsquares = int(sys.argv[1])

xlim=([-numpy.pi,numpy.pi])
ylim=([-numpy.pi,numpy.pi])

print xlim[0]
print xlim[1]
print ylim

midpts = numpy.zeros((numgridsquares,numgridsquares,2))
boundaries = numpy.zeros((numgridsquares, numgridsquares, 4))

print midpts[1,1]
print "TREV"
xdiff = (xlim[1]-xlim[0])/numgridsquares

print xdiff

for i in range(0,numgridsquares):
    tmpx = xlim[0]+xdiff*(i+0.5)
    for j in range(0,numgridsquares):
        tmpy = ylim[0]+xdiff*(j+0.5)
        print "tmpx ",tmpx,", tmpy ",tmpy 
        midpts[i][j][0] = tmpx
        midpts[i][j][1] = tmpy
        boundaries[i][j][0] = xlim[0]+xdiff*i
        boundaries[i][j][1] = ylim[0]+xdiff*j
        boundaries[i][j][2] = xlim[0]+xdiff*(i+1)
        boundaries[i][j][3] = ylim[0]+xdiff*(j+1)

print midpts[0][2][0]," ",midpts[0][2][1]
print "A!)"

physarray=[]
ioarray=[]
proparray=[]
forces=[]


for i in range(0,numgridsquares*numgridsquares):
    physarray.append(Physical())
    ioarray.append(IO())
    forces.append(Forces())
    ioarray[i].readPDBPos(physarray[i],"data/mydata/minC7eq.pdb")
    ioarray[i].readPSF(physarray[i],"data/mydata/alan_mineq.psf")
    ioarray[i].readPAR(physarray[i],"data/mydata/par_all27_prot_lipid.inp")
    physarray[i].bc = "Vacuum"
    physarray[i].temperature=300
    physarray[i].exclude="scaled1-4"
    physarray[i].seed = 1234
    proparray.append(Propagator(physarray[i],forces[i],ioarray[i]))


ff = forces[0].makeForceField(physarray[0])
print "AA"
#ff.bondedForces("badihh")
ff.nonbondedForces("lc")

kappa=20.0
PHI=11
PSI=18

#ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
#                                 'dihedralnum':[PHI-1, PSI-1],
#                                 'angle':[midpts[0][0][0],midpts[0][0][1]]}


print "G"
#stringgraph=ioarray[0].newGraph('Phi','Psi')
bins = {}
for ii in range(0,numgridsquares):
    for jj in range(0,numgridsquares):

        ff.bondedForces("badihh")
        ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                         'dihedralnum':[PHI-1, PSI-1],
                                         'angle':[midpts[0][0][0],midpts[0][0][1]]}
        forces[0].build()
        #if (ii !=0):
        #    setPhiPsi(midpts[ii][jj][0],midpts[ii][jj][1],ff)
        #elif (jj != 0) :
        #    setPhiPsi(midpts[ii][jj][0],midpts[ii][jj][1],ff)

        #raw_input()
        print "H"
        pts = ii*numgridsquares + jj
        print pts
        ef = 'data/movetomidpt.energies.'+str(ii)+'.'+str(jj)
        df = 'data/movetomidpt.'+str(ii)+'.'+str(jj)+'.dcd'
        print ef," ",df
        #ioarray[pts].files={'energies':(ef,100),'dcdtraj':(df,10)}

        for i in range(0,10):
            proparray[pts].propagate(scheme="LangevinImpulse",steps=1000,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
            #ioarray[pts].plotVector(proparray[pts],stringgraph,[[physarray[pts].angle(PHI),physarray[pts].angle(PSI)]],rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi]) 
        

        ff.bondedForces('badihhhhhhhh')
        ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa, kappa, kappa, kappa, kappa, kappa, kappa],
                                         'dihedralnum':[PHI-1, PSI-1, PHI-1, PSI-1, PHI-1, PSI-1, PHI-1, PSI-1],
                                         'angle':[boundaries[ii][jj][0],boundaries[ii][jj][1],boundaries[ii][jj][0],boundaries[ii][jj][3],
                                                  boundaries[ii][jj][2],boundaries[ii][jj][1],boundaries[ii][jj][2],boundaries[ii][jj][3]]}
        forces[0].build()
        for i in range(0, 200000):
            proparray[pts].propagate(scheme="LangevinImpulse",steps=1,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
            #ioarray[pts].plotVector(proparray[pts],stringgraph,[[physarray[pts].angle(PHI),physarray[pts].angle(PSI)]],rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
            if (not (bins.has_key(str(physarray[pts].angle(PHI))+" "+str(physarray[pts].angle(PSI))))):
               bins[str(physarray[pts].angle(PHI))+" "+str(physarray[pts].angle(PSI))] = 1
            else:
               bins[str(physarray[pts].angle(PHI))+" "+str(physarray[pts].angle(PSI))] += 1

bins = [(v, k) for k, v in bins.items()]
bins.sort()
bins.reverse()
bins = [(k, v) for v, k in bins]

for aa in range(0, len(bins)):
    print bins[aa][0], " ", bins[aa][1]

#sys.exit(0)

        #.files={'energies':('contour.energies.1',10),'dcdtraj':('contour.1.dcd',10)}
#myIO.screen=1
#stringgraph=myIO.newGraph('Phi', 'Psi');
#totalsteps=10000
#z=[]
#for i in range(0,20):
    #myProp.propagate(scheme="LangevinImpulse",steps=100,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
    #z.append([myPhys.angle(PHI),myPhys.angle(PSI)])
    #myIO.plotVector(myProp,stringgraph,z,rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
    #z.pop()
