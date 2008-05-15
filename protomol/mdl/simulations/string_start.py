from Physical import *
from Forces import *
from Propagator import *
from IO import *
from ForceField import *

import math
import numpy
import string

####Call this program like this :
####python string_start.py <debuglevel (1 is preferred)> <common vectors> <init string outfile> <eigvecfile>

commandlineargs = len(sys.argv)
debug = 0

if (commandlineargs > 1):
    debug = int(sys.argv[1])

#no of common vectors
if (commandlineargs > 2):
    cv = int(sys.argv[2])
else:
    print 'Using all vectors from eigenvector file'
    cv = -1

#3rd command line argument
#output filename for the initial string
if (commandlineargs > 3):
    initialstringoutputfilename = sys.argv[3]

if (commandlineargs > 4):
    eigvecfilename = sys.argv[4]
else :
    eigvecfilename = "data/eigVmC7mwtxt"

io = IO()

numpoints = 10

#If In have r points on the string, I need r physical objects.

physarray = []
forcearray = []
proparray = []


#Initialization routine

def string_init():
    for i in range(0,numpoints):
        physarray.append(Physical())
        forcearray.append(Forces())
        if (i > numpoints / 2):
            io.readPDBPos(physarray[i],"data/armin_a3.pdb")
            #io.readPDBPos(physarray[i],"data/minC7eq_a3.pdb")
        else:
            io.readPDBPos(physarray[i],"data/minC7eq_a3.pdb")
            #io.readPDBPos(physarray[i],"data/armin_a3.pdb")
        io.readPSF(physarray[i], "data/alan_mineq.psf")
        io.readPAR(physarray[i], "data/par_all27_prot_lipid.inp")

        #Read the eigenvector file
        io.readTextEigenvectors(physarray[i], eigvecfilename)
        physarray[i].bc = "Vacuum"
        physarray[i].temperature = 300
        physarray[i].exclude = "scaled1-4"
        physarray[i].seed = 1234
        if(debug==10):
            print (physarray[i].angle(11)*180)/numpy.pi
            print (physarray[i].angle(18)*180)/numpy.pi
        proparray.append(Propagator(physarray[i],forcearray[i],io))

def sqrtmass(rd,cd,md):
    mmat = numpy.ndarray([rd,cd])
    #print mmat.size
    for i in range(0,rd):
        for j in range(0,cd):
            if (i==j):
                mmat[i,j] = numpy.sqrt(md[i])
            else:
                mmat[i,j] = 0
    return mmat

def matmult6(m,v):
    w = m.__mul__(v)
    return w

def convert_to_numpy_matrix(ma,colsize):
    numcols = ma.size/colsize
    mv = numpy.zeros((colsize,numcols))
    for j in range(0,numcols):
        for i in range(0,colsize):
            mv[i,j] = ma[j*colsize+i]

    return mv
            
def getccols(ccol,csu,j):
    n = csu.size
    #print 'n = ',n
    for i in range(0,n):
        ccol[i,j] = ccol[i,j-1]+csu[i,0]

def mvmul(m,v):
    r1 = m.size/m[0].size
    r2 = m[0].size
    #print 'r1 ',r1
    #print 'r2 ',r2
    mv = numpy.zeros((r1,1))
    for i in range(0,r1):
        sum = 0
        for j in range(0,r2):
            sum = sum + m[i,j]*v[j]
        mv[i,0]=sum
    return mv

def addv(x,y,v):
    #v = x.copy()
    r1 = y.size
    #print 'y size ',r1
    for i in range(0,r1):
        v[i] = v[i]+y[i] 

#phys stands for a physical object and
#x contains the positions I need to update
def copy_phys_positions(phys,x):
   for ii in range(0,x.size):
       phys.positions[ii]=x[ii]
    

        

# Main code block
string_init()


######## Calculate c_i for each point on the string

##Get the Q_S matrix
xa = physarray[0].positions
xb = physarray[numpoints-1].positions
#print 'Position a'
#print xa
#print 'Position b'
#print xb

#myeigvecs = physarray[0].myEig;

Q_S = io.getTextEigenvectorsData(physarray[0])

###%%%%%%%%%%%%% Debug statements %%%%%%%%%%%###############
#print Q_S.size
#print Q_S[0]
#print Q_S[65]
#print Q_S[66]
#print Q_S[131]
#print xa.size
###%%%%%%%%%%%%% EndDebug statements %%%%%%%%%##############

######### I need to take only those column vectors which are common ###############

Q_S_M = convert_to_numpy_matrix(Q_S,xa.size)
#print 'Final eigenvector matrix :'
#print Q_S_M

###%%%%%%%%%%%%% Debug statements %%%%%%%%%%%###############
#if(debug >5):
    #print Q_S_M[0,0]
    #print Q_S_M[65,0]
    #print Q_S_M[0,1]
    #print Q_S_M[65,1]
    #print Q_S_M[0,2]
    #print Q_S_M[65,2]
    #print Q_S_M.size
    #print Q_S_M[:,3]
    #print Q_S_M
###%%%%%%%%%%%%% EndDebug statements %%%%%%%%%##############

if (cv == -1):
    cv = Q_S_M.size/Q_S_M[0].size

######### I need to take only those column vectors which are common ###############
Q_S_M = Q_S_M[:,0:cv]
#print 'Final eigenvector matrix :'
#print Q_S_M
#if (debug>=1):
    #print Q_S_M.size
    #print Q_S_M[0].size

##Get the matrix M^{1/2}
sqmH = sqrtmass(len(physarray[0].masses),len(physarray[0].masses),physarray[0].masses)
xdf = xb-xa


####### SC : Next two lines generate the c vector from 0 to 1 along the string #######
#temp = matmult6(numpy.matrix(sqmH),numpy.matrix(xdf).transpose())
#c_s = matmult6(numpy.matrix(Q_S_M).transpose(),temp)
if(debug>5):
    print 'xdf '
    print xdf

if(debug>6):
    print 'sqmH matrix '
    print sqmH

#temp_m = numpy.matrix(numpy.matrix(sqmH) * numpy.matrix(xdf).transpose())
temp_m = mvmul(sqmH,xdf)
if(debug>5):
    print 'tempm matrix '
    print temp_m
#c_s = numpy.matrix(Q_S_M.transpose() * temp_m)
c_s = mvmul(Q_S_M.transpose(),temp_m)
if(debug>5):
    print 'c_s matrix '
    print c_s

####### SC: Diff in c value between contiguous points on the string #################
c_s_u = c_s/(numpoints-1)

if(debug>5):
    print 'c_s_u size ',c_s_u.size
    print c_s_u
####### SC: Initializing the 2D array for storing c values at each point on the string #####
cmodevalarray = numpy.zeros((c_s_u.size,numpoints))
#print cmodevalarray[:,0]
#getccols(cmodevalarray,c_s_u,1)
#print cmodevalarray[:,1]

####### SC: fill in cmodevalarray matrix ############
for i in range(1,numpoints):
    getccols(cmodevalarray,c_s_u,i)

#print 'Full cmodevalarray output : '
#print cmodevalarray
#print 'cmodevalarray last column : '
#print cmodevalarray[:,numpoints-1]

if (debug>5):
    print 'cmodevalarray 1' 
    print cmodevalarray[:,1]
    print 'cmodevalarray 2' 
    print cmodevalarray[:,2]

invsqmH = numpy.linalg.inv(sqmH)

#print 'invsqmH : '
#print invsqmH

for i in range(0,numpoints):
    xx = mvmul(Q_S_M,cmodevalarray[:,i])
    xf = mvmul(invsqmH,xx)
    xp = xa.copy()
    addv(xa,xf,xp)
    if (debug>5):
        print 'Old positions '
        print physarray[i].positions
        print 'New positions '
    #physarray[i].positions = xp.copy()
    #physarray[i].__setattr__('positions',xp.copy())
    copy_phys_positions(physarray[i],xp.copy())
    fn = 'spt'+str(i)+'.pdb'
    if (debug>5):
        print fn
    #print xp.size
    #print xp[0].size
    io.writePDBPos(physarray[i],fn)
    if (debug>5):
        print xp

    ###### SC : Just checking
    xdiff = xp - xa;
    cdiff = mvmul(Q_S_M.transpose(),mvmul(sqmH,xdiff))
    #print 'C vals after coord move : '
    #print cdiff



####### SC: Initialize propagator for 0 steps in order to get reference to it ##########
for ii in range(1,numpoints-1):
    ###### SC: Define force fields ##################
    ff1 = forcearray[ii].makeForceField(physarray[ii])
    ff1.bondedForces("badi")
    ff1.nonbondedForces("le")

    ff2 = forcearray[ii].makeForceField(physarray[ii])
    ff2.bondedForces("badi")
    ff2.nonbondedForces("le")
    print 'Starting equilibration at point '+str(ii)

    io.screen = 1
    proparray[ii].propagate(scheme=["NormModeResInt", "NormModeResMin"],
       steps=1000, cyclelength=1, dt=1.0, forcefield=[ff1, ff2],
       params={'NormModeResInt':{'fixmodes':52,'gamma':91,'fdof':0, 'kval':100},
               'NormModeResMin':{'avModeMass':30}})


#for ii in range(1,numpoints-1):
#    proparray[ii].myPropagator.setDataRefC(cmodevalarray[:,ii])




if(debug == 1):
    if (commandlineargs > 3):
        initialstringoutputfile = open(initialstringoutputfilename,'w')
        
    for i in range(0,numpoints):
        phi = (physarray[i].angle(11)*180)/numpy.pi
        psi = (physarray[i].angle(18)*180)/numpy.pi
        s = str(phi)+' '+str(psi)+'\n'

        if (commandlineargs > 3):
            initialstringoutputfile.write(s)
        #else :
            #print 'Phi Psi ',i
            #print phi,' ',psi

    if (commandlineargs > 3):
        initialstringoutputfile.close()



