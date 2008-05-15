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

print debug

#no of common vectors
cv = -1
if (commandlineargs > 2):
    cv = int(sys.argv[2])

print cv

#3rd command line argument
#output filename for the initial string
if (commandlineargs > 3):
    initialstringoutputfilename = sys.argv[3]
    # delete this file if its already there
    cmd = 'rm -f '+initialstringoutputfilename
    os.system(cmd)

if (commandlineargs > 4):
    eigvecfilename = sys.argv[4]


if (commandlineargs > 5):
	numpoints = int(sys.argv[5])
else :
	numpoints = 25

# eigenvalue file name
if (commandlineargs > 6):
    eigvalfilename = sys.argv[6]



#io = IO()


#If In have r points on the string, I need r physical objects.

physarray = []
forcearray = []
proparray = []
ioArray = []

for iii in range(0,numpoints):
    ioArray.append(IO())

#kval=10

gamma=500.0

totalsteps = 1000

#Initialization routine

def string_init():
    for i in range(0,numpoints):
        #print "####################### Initializing point "+str(i)+" ########################"
        #print "##############################################################################"
        #print "Before physical constructor"
        physarray.append(Physical())
        #print "After physical constructor"
        forcearray.append(Forces())

        ##### Every time I update an element of physical array, I might need to rebuild
        ##### the physical pbject.

        if (i > numpoints / 2):
            ioArray[i].readPDBPos(physarray[i],"data/armin_a3.pdb")
            #io.readPDBPos(physarray[i],"data/minC7eq_a3.pdb")
        else:
            ioArray[i].readPDBPos(physarray[i],"data/minC7eq_a3.pdb")
            #io.readPDBPos(physarray[i],"data/armin_a3.pdb")
        ioArray[i].readPSF(physarray[i], "data/alan_mineq.psf")
        ioArray[i].readPAR(physarray[i], "data/par_all27_prot_lipid.inp")

        #Read the eigenvector file
        ioArray[i].readTextEigenvectors(physarray[i], eigvecfilename)
        physarray[i].bc = "Vacuum"
        physarray[i].temperature = 300
        physarray[i].exclude = "scaled1-4"
        physarray[i].seed = 1234
        if(debug==10):
            print (physarray[i].angle(11)*180)/numpy.pi
            print (physarray[i].angle(18)*180)/numpy.pi

        ef = 'energies.out.'+str(i)
        ioArray[i].files = {'energies':(ef,1)}
        proparray.append(Propagator(physarray[i],forcearray[i],ioArray[i]))

os.system('rm -f energies.out.*')

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
            #print i," ",j
            sum = sum + m[i,j]*v[j]
        mv[i,0]=sum
    return mv

def mmul(m1,m2):
    m1_1 = m1.size/m1[0].size
    m1_2 = m1[0].size
    #print 'm1_1 ',m1_1
    #print 'm1_2 ',m1_2
    
    m2_1 = m2.size/m2[0].size
    m2_2 = m2[0].size
    #print 'm2_1 ',m2_1
    #print 'm2_2 ',m2_2
   
    mv = numpy.zeros((m1_1,m2_2)) 
    for i in range(0,m1_1):
        for j in range(0,m2_2):
            for k in range(0,m1_2):
                mv[i,j] = mv[i,j] + m1[i,k]*m2[k,j]

    return mv
    

def addv(y,v):
    r1 = y.size
    #print 'y size ',r1
    for i in range(0,r1):
        v[i] = v[i]+y[i] 

#phys stands for a physical object and
#x contains the positions I need to update
def copy_phys_positions(phys,x):
    for ii in range(0,x.size):
        phys.positions[ii]=x[ii]


#
# Updates new string in cmarray
#
def updateString(myPropagatorClass, cmarray, ii, vsteps, dt):
    newC = myPropagatorClass.myPropagator.getData()
    #kvals = myPropagatorClass.myPropagator.getKvals()
    #print newC
    tp = (newC-cmarray[:,ii])
    #print "cdiff size .... ",(newC-cmarray[:,ii]).size
    tempArrayii = eigvals.transpose()*(newC-cmarray[:,ii])
    print "tempArrayii size ",tempArrayii.size
    #print tempArrayii
    tempArrayii = mvmul(numpy.transpose(Mtd),tempArrayii.transpose())
    #tempArrayii = mmul(numpy.transpose(Mtd),tempArrayii.transpose())
    #print "iith cmodevalarray size ",cmarray[:,ii].size
    xxx = (tempArrayii*vsteps*dt)/gamma
    #print xxx
    #print tempArrayii
    addv(cmarray[:,ii],(-1)*xxx)


def EigvalReader(eigvalfilename):
    myeigfile = open(eigvalfilename,"r")
    evals = numpy.zeros((cv,1))
    i= 0
    for line in myeigfile:
        evals[i] = float(line)
        i = i+1
    myeigfile.close()
    for ii in range(0,evals.size):
        evals[ii] = 5*evals[ii]
        if (evals[ii] < 1):
            evals[ii] = 1

    return evals


# Calculates total length until point p on the string
def CalcNorm(V,p):
    result = 0
    if(p==0):
        return result
    for i in range(1,p):
        #print i
        #print V[:,i]
        #print V[:,i-1] 
        #result += numpy.linalg.norm(addv(V[:,i],(-1)*V[:,i-1]))
        result += numpy.linalg.norm(V[:,i]-V[:,i-1])
    return result
    #L = numpy.linalg.norm(V) 

def lval(V,p):
    #R = len(V)
    R = numpoints
    #print "p c ",R
    return p*CalcNorm(V,R)/(R-1)

def qval(V,p):
    tmpQ = 1
    #print "p b ",p
    myl = lval(V,p)
    #print "myl ",myl
    while (not (CalcNorm(V,tmpQ-1) < myl and myl<=CalcNorm(V,tmpQ))):
        tmpQ += 1
        #print "tmpQ ",tmpQ
    return tmpQ

def SmoothingAndReparameterization(cm):
    S = []
    S.append(cm[:,0])
    #print cm.size
    #print cm[:,0].size
    for p in range(1,numpoints-1):
        #print "p a ",p
        myQ = qval(cm,p)
        newCatp = cm[:,myQ] + (lval(cm,p)-CalcNorm(cm,myQ-1))*(cm[:,myQ]-cm[:,myQ-1])/numpy.linalg.norm(cm[:,myQ]-cm[:,myQ-1])
        S.append(newCatp)
    S.append(cm[:,numpoints-1])
    return S
   


#
# My new propagate fuction which overrides the default propagate function.
#
def simplePropagateFuction(myPropagatorClass, cmarray, ii, vsteps, scheme="Leapfrog", steps=0, cyclelength=-1, dt=0.1, forcefield=[], params={}):
    myPropagatorClass.myTimestep = dt
    chain = ()
    if (cyclelength != -1):  # MTS
        if (str(type(cyclelength))[7:11] == 'list'): # LIST, MANY LEVELS
            levels = len(cyclelength) + 1
            outertime = cyclelength[0]
        else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
            levels = 2
            outertime = cyclelength

        if (str(type(scheme))[7:11] == 'list'): # LIST, MANY LEVELS
            outerscheme = scheme[0]
            #print outerscheme
        else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
            outerscheme = scheme


        # THE NUMBER OF FORCEFIELDS PROVIDED MUST EQUAL THE NUMBER
        # OF PROPAGATION LEVELS
        if (len(forcefield) != levels):
            print "[MDL] Error in propagate(): ", levels, " levels of propagation with ", len(forcefield), " force fields."
        outerforcefield = forcefield[0]

        if (str(type(scheme))[7:11] != 'list'):
            chain += (params,)
        else:
            if (params.has_key(outerscheme)):
                chain += (params[outerscheme],)
            else:
                chain += ({},)

        for i in range(1, levels):
            if (str(type(scheme))[7:11] == 'list' and i < len(scheme)):
                chain += (scheme[i],)
            if (str(type(cyclelength))[7:11] == 'list' and i < len(cyclelength)):
                chain += (cyclelength[i],)
            else:
                chain += (dt,)
            chain += (forcefield[i],)
            if params.has_key(scheme[i]):
                chain += (params[scheme[i]],)
            else:
                chain += ({},)

        #print chain


    else: #STS
        outertime = dt
        outerscheme = scheme
        outerforcefield = forcefield
        chain += (params,)

    if (myPropagatorClass.forces.dirty()):
        myPropagatorClass.forces.build()

    if (myPropagatorClass.io.dirty):
        myPropagatorClass.io.build()

    #print outerscheme
    if (propFactory.getType(outerscheme) == "method"):
        # Calculate the forces, store them in force.
        print 'Calculate initial force for outer scheme '
        outerforcefield.calculateForces(myPropagatorClass.phys, myPropagatorClass.forces)
        print 'Initialize : Updating center of mass and momenta '
        myPropagatorClass.phys.updateCOM_Momenta()
        print 'Initialize : Run IO module to generate initial output '
        myPropagatorClass.io.run(myPropagatorClass.phys, myPropagatorClass.forces, 0, outertime)
        myPropagatorClass.io.myProp = myPropagatorClass
        for ii in range(1, steps+1):
            propFactory.create(outerscheme, myPropagatorClass.phys, myPropagatorClass.forces, myPropagatorClass.io, 1, outertime*Constants.invTimeFactor(), outerforcefield, *chain)
            #if (ii % 1000 == 1):
            #   print (ii-1)/1000, " ps"
            #self.phys.time = ii*outertime*Constants.invTimeFactor()
            #self.io.run(self.phys, self.forces, ii, outertime)
            #self.phys.updateCOM_Momenta()
    else:
        #print "Propagator is an object"
        #sys.exit(0)
        #return
        #print chain
        setPropagator(myPropagatorClass, myPropagatorClass.phys, myPropagatorClass.forces, propFactory.applyModifiers(propFactory.create(outerscheme, outertime, outerforcefield, *chain), outerscheme))
        ## Updating the eigenvalues/rayleigh quotients
        myPropagatorClass.myPropagator.SetEigvals(eigvals)

        ## Estimate k values for each mode
        #myPropagatorClass.myPropagator.SetKvals()
        xa = physarray[0].positions
        myPropagatorClass.myPropagator.SetInitialPos(xa.size,xa)
        myPropagatorClass.myPropagator.copyToVector3DBlock()

        shake = False
        if (params.has_key('shake') and params['shake'] == 'on'):
            shake = True
            shakeMod = myPropagatorClass.myPropagator.createShakeModifier(0.00001, 30)
            myPropagatorClass.myPropagator.adoptPostDriftOrNextModifier(shakeMod)
        rattle = False
        if (params.has_key('rattle') and params['rattle'] == 'on'):
            rattle = True
            rattleMod = myPropagatorClass.myPropagator.createRattleModifier(0.02, 30)
            myPropagatorClass.myPropagator.adoptPostDriftOrNextModifier(rattleMod)

        # call setDataRefC
        myPropagatorClass.myPropagator.setDataRefC(cmarray[:,ii])
        executePropagator(myPropagatorClass, myPropagatorClass.phys, myPropagatorClass.forces, myPropagatorClass.io, vsteps)

        updateString(myPropagatorClass,cmarray,ii,vsteps,dt)

        cmarray = SmoothingAndReparameterization(cmarray)



def GetQSMatrix(c_v):
    xa = physarray[0].positions
    Q_S = ioArray[0].getTextEigenvectorsData(physarray[0])
    Q_S_M = convert_to_numpy_matrix(Q_S,xa.size)
    if (c_v == -1):
        c_v = Q_S_M.size/Q_S_M[0].size

    ######### I need to take only those column vectors which are common ###############
    Q_S_M = Q_S_M[:,0:c_v]
    return Q_S_M

def EndToEndModeDiff(Q_S_M,sqmH):
    xa = physarray[0].positions
    xb = physarray[numpoints-1].positions
    xdf = xb-xa


    ####### SC : debugging... #######
    if(debug>5):
        print 'xdf '
        print xdf

    if(debug>6):
        print 'sqmH matrix '
        print sqmH

    #### end debug ##################

    #### temp_m = M^{1/2}(x-x0)
    temp_m = mvmul(sqmH,xdf)

    ####### SC : debugging ... #######
    if(debug>5):
        print 'tempm matrix '
        print temp_m
    ####### end debug ################


    #### c_s = Q_S^{T}M^{1/2}(x-x0)
    c_s = mvmul(Q_S_M.transpose(),temp_m)


    ####### SC : debugging ... #######
    if(debug>5):
        print 'c_s matrix '
        print c_s
    ####### end debug ################
    return c_s

def Mtensor(Q_S_M,sqmH,invsqmH):
	#####Calculating M tensor
    Mtd = mmul(numpy.transpose(Q_S_M),sqmH)
    invM = mmul(invsqmH,invsqmH)
    MM = mmul(invM,numpy.transpose(Mtd))
    Mtd = mmul(Mtd,MM)
    return Mtd

   
def Equilibrate_string(): 
    for ii in range(1,numpoints):
        ff1 = forcearray[ii].makeForceField(physarray[ii])
        ff1.bondedForces("badi")
        ff1.nonbondedForces("le")

        ff2 = forcearray[ii].makeForceField(physarray[ii])
        ff2.bondedForces("badi")
        ff2.nonbondedForces("le")
        print 'Starting equilibration at point '+str(ii)

        #print "physarray "+str(ii)+" velocities "
        #print physarray[ii].velocities
        ioArray[ii].screen = 1
        proparray[ii].propagate(scheme=["NormModeResInt", "NormModeResMin"],
            steps=25, cyclelength=1, dt=1.0, forcefield=[ff1, ff2],
            params={'NormModeResInt':{'fixmodes':66-cv,'gamma':91,'fdof':0},
            'NormModeResMin':{'avModeMass':3.0}})


def OutputStringPos():
    initialstringoutputfile = open(initialstringoutputfilename,'a')

    z = []
    for i in range(0,numpoints):
        phi = (physarray[i].angle(11)*180)/numpy.pi
        psi = (physarray[i].angle(18)*180)/numpy.pi
        z.append([phi,psi])
        s = str(phi)+' '+str(psi)+'\n'
        initialstringoutputfile.write(s)

    initialstringoutputfile.close()
    return z



################### Main code block ################################################

string_init()


##Get the Q_S matrix
Q_S_M = GetQSMatrix(cv)


##Get the matrix M^{1/2}
sqmH = sqrtmass(len(physarray[0].masses),len(physarray[0].masses),physarray[0].masses)

c_s = EndToEndModeDiff(Q_S_M,sqmH)

####### SC: Diff in c value between neighboring points on the string #################
c_s_u = c_s/(numpoints-1)

####### SC : debugging ... #######
if(debug>5):
    print 'c_s_u size ',c_s_u.size
    print c_s_u
####### end debug ################


####### SC: Initializing the 2D array for storing c values at each point on the string #####
cmodevalarray = numpy.zeros((c_s_u.size,numpoints))

####### SC: fill in cmodevalarray matrix ############
for i in range(1,numpoints):
    getccols(cmodevalarray,c_s_u,i)

print cmodevalarray.size

print cmodevalarray[0].size
print cmodevalarray[:,0].size

#sys.exit(0)

if (debug>5):
    print 'cmodevalarray 1' 
    print cmodevalarray[:,1]
    print 'cmodevalarray 2' 
    print cmodevalarray[:,2]

invsqmH = numpy.linalg.inv(sqmH)

Mtd = Mtensor(Q_S_M,sqmH,invsqmH)

#sys.exit(0)

xa = physarray[0].positions

for i in range(0,numpoints-1):
    xx = mvmul(Q_S_M,cmodevalarray[:,i])
    xf = mvmul(invsqmH,xx)
    xp = xa.copy()
    addv(xf,xp)
    if (debug>5):
        print 'Old positions '
        print physarray[i].positions
        print 'New positions '
    copy_phys_positions(physarray[i],xp.copy())
    fn = 'spt'+str(i)+'.pdb'
    if (debug>5):
        print fn
    ioArray[i].writePDBPos(physarray[i],fn)
    if (debug>5):
        print xp

    ###### SC : Just checking
    xdiff = xp - xa;
    cdiff = mvmul(Q_S_M.transpose(),mvmul(sqmH,xdiff))
    #print 'C vals after coord move : '
    #print cdiff

z=OutputStringPos()

Equilibrate_string()

z=OutputStringPos()

eigvals = EigvalReader(eigvalfilename)

#myIO = IO()
stringgraph=ioArray[0].newGraph('Phi', 'Psi')

####### SC: Initialize propagator and minimize for few steps ##########
for k in range(0,totalsteps):
    for ii in range(1,numpoints-1):
        ###### SC: Define force fields ##################
        ff1 = forcearray[ii].makeForceField(physarray[ii])
        ff1.bondedForces("badi")
        ff1.nonbondedForces("le")

        ff2 = forcearray[ii].makeForceField(physarray[ii])
        ff2.bondedForces("badi")
        ff2.nonbondedForces("le")
        #print 'Starting equilibration at point '+str(ii)

        ioArray[ii].screen = 1
        simplePropagateFuction(proparray[ii],cmodevalarray,ii,1,scheme=["NormModeSampling", "NormModeSamplingMin"],
            steps=1, cyclelength=1, dt=1.0, forcefield=[ff1, ff2],
            params={'NormModeSampling':{'fixmodes':66-cv,'gamma':91,'fdof':0,'eqlb':0},
               'NormModeSamplingMin':{'avModeMass':3.0}})

    z = OutputStringPos()
    print z
    #ioArray[0].plotVector(proparray[0],stringgraph,z,rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])
    ioArray[0].plotVector(proparray[0],stringgraph,z,rangex=[-180.0, 180.0], rangey=[-180.0, 180.0])
