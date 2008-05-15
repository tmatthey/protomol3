from __future__ import generators
from Vector3DBlock import *
import math
import topologyutilities
import numpy
import sys

class SplineInterpolator:
	"""
	Cubic spline interpolator.
	"""
	def __init__(self):
		#print "Spline interpolator initialized"
		self.diag = []     #: Middle diagonal elements
		self.u_diag = []   #: Upper diagonal elements
		self.l_diag = []   #: Lower diagonal elements
		self.R = []        #: Right hand side
		self.L = []        #: Lower bidiagonal matrix
		self.U = []        #: Upper bidiagonal matrix
		self.xx = []       #: Holds solution

	###############################################################
	#
	#   diag stands for the diagonal elements of tridiagonal matrix
	#   u_diag is the upper diagonal element (f_i in van loan book)
	#   l_diag : lower diagonal elements.
	#   In order to keep the size of diag, u_diag and l_diag same,
	#   I add 0 at the start or end (as necessary).
	#
	###############################################################

	def SetTridiagonalSystem(self,x,y,num):
		"""
		Set up the tridiagonal system

		@type x: list
		@param x: x values

		@type y: list
		@param y: y values

		@type num: int
		@param num: Number of points
		"""
		for i in range(0,num-2):
        	#print i
			dxi = x[i+1]-x[i]
			dxi1 = x[i+2]-x[i+1]
			dyi = (y[i+1]-y[i])/dxi
			dyi1 = (y[i+2]-y[i+1])/dxi1

			if i == 0 :
				self.diag.append( 2*dxi+dxi1+((dxi*dxi)/dxi1))
				self.u_diag.append(dxi+((dxi*dxi)/dxi1))
				self.l_diag.append(0)
				self.R.append((dxi1*dyi+3*dxi*dyi1 + 2*dyi1*((dxi*dxi)/dxi1)))
			elif i == (num-3):
				self.u_diag.append(0)
				self.l_diag.append(dxi1+((dxi1*dxi1)/dxi))
				self.diag.append(2*dxi1+dxi+((dxi1*dxi1)/dxi))
				self.R.append(3*dxi1*dyi+dxi*dyi1 +2*dyi*((dxi*dxi)/dxi1))
			else:
				self.diag.append(2*(dxi+dxi1))
				self.l_diag.append(dxi1)
				self.u_diag.append(dxi)
				self.R.append(3*(dxi1*dyi+dxi*dyi1))

		##################################################
		#  Just for debugging
		##################################################
		#print self.diag
		#print self.u_diag
		#print self.l_diag

	###################################################################
	#
	#  Construct two bi-diagonal matrices L and U.
	#
	###################################################################

	def FindLU(self,num):
		"""
		Construct two bi-diagonal matrices L and U.

		@type num: int
		@param num: Number of points
		"""
		for i in range(0,num):
			if i == 0:
				self.U.append(self.diag[i])
				self.L.append(0)
			else:
				self.L.append(self.l_diag[i]/self.U[i-1])
				l_temp = self.l_diag[i]/self.U[i-1]
				self.U.append(self.diag[i] - l_temp*self.u_diag[i])

	####################################################################
	#
	# LBiDiSol : Lower Bi-diagonal solver.
	#            Solves lower bi-diagonal system Lx = b
	#
	####################################################################

	def LBiDiSol(self,num):
		"""
		Lower Bi-diagonal solver.
		Solves lower bi-diagonal system Lx = b

		@type num: int
		@param num: Number of points
		"""
		self.xx.append(self.R[0])
		for i in range(1,num):
			self.xx.append(self.R[i]-self.L[i]*self.xx[i-1])

	####################################################################
	#
	#  UBiDiSol : Upper bi-diagonal solver.
	#
	####################################################################
	
	def UBiDiSol(self,num):
		"""
		Upper Bi-diagonal solver.

		@type num: int
		@param num: Number of points

		@rtype: list
		@return: Solution to uper bi-diagonal system
		"""
		S = []
		for i in range(0,num+2):
			S.append(0)
		S[num-1] = self.R[num-1]/self.U[num-1]
		for i in range(num-2,-1,-1):
			S[i] = (self.R[i] - self.u_diag[i]*S[i+1])/self.U[i]

		return S

	def locate(self,x,k):
		"""
		Locate where in the x points the passed value falls

		@type x: list
		@param x: List of x values

		@type k: float
		@param k: Value

		@rtype: int
		@return Index into x values
		"""
                #print k
                #print x
		xs = x.__len__()
		for i in range(1,xs):
			if k <= x[i]:
				return i-1

	#####################################################################
	#
	# Evaluate approximated function at values nx.
	#
	#####################################################################

	def EvalSpline(self,x,y,S,nx):
		"""
		Evaluate approximated function at values nx

		@type x: list
		@param x: X values

		@type y: list
		@param y: Y values

		@type S: list
		@param S: Appropximated function

		@type nx: list
		@param nx: Values at which to evaluate

		@rtype: list
		@return: Approximated function values
		"""
		nxs = nx.__len__()
		#print "Size of array for which function to be evaluated "
		#print nxs
	
		nf = []
		for i in range(0,nxs):
			# Locates the interval where nx[i] lies 
			loc = self.locate(x,nx[i])
			if loc == x.__len__()-1 :
				nf.append(x[x.__len__()-1])
				break
			#print "index "+str(loc)+" and i "+str(i)
			dxi = x[loc+1]-x[loc]
			yprime = (y[loc+1]-y[loc])/dxi
			zmxi = nx[i]-x[loc]
			fval = y[loc]+S[loc]*zmxi+((yprime-S[loc])/(dxi))*(zmxi*zmxi)+((S[loc]+S[loc+1]-2*yprime)/(dxi*dxi))*(zmxi*zmxi)*(nx[i]-x[loc+1]);
			nf.append(fval)

		return nf

	def Spline(self,x,y):
		"""
		Highest level routine for setting up member matrices.

		@type x: list
		@param x: X values

		@type y: list
		@param y: Y values

		@rtype: list
		@return: Interpolated values
		"""
		num = x.__len__()
		self.SetTridiagonalSystem(x,y,num)
		self.FindLU(num-2)
		self.LBiDiSol(num-2)
		S = self.UBiDiSol(num-2)
		return S

	def SmoothingAndReparameterization(self,x,y,xx,yy,num):
		"""
		Reparameterize by arc length

		@type x: list
		@param x: X values

		@type y: list
		@param y: Y values

		@type xx: list
		@param xx: Linear space

		@type yy: list
		@param yy: Corresponding spline values

		@type num: int
		@param num: Number of (x,y) points

		@rtype: list
		@return: (x,y) points after reparameterization
		"""
		len = []
		sx = xx.__len__()
		for i in range(0,sx-1):
			len.append(math.sqrt((xx[i]-xx[i+1])**2 + (yy[i]-yy[i+1])**2))
			if (i > 0):
				len[i] += len[i-1]

		stp = len[sx-2]/(num-1)
		stpcmp=0
		xo = []
		yo = []
		xo.append(x[0])
		yo.append(y[0])
		nextj=0
		for ln in range(1, num-1):
			stpcmp += stp
			for j in range(nextj, sx-1):
				if (len[j] > stpcmp):
					xo.append(xx[j])
					yo.append(yy[j])
					nextj = j
					break
		xo.append(x[num-1])
		yo.append(y[num-1])
		retval = []
		for i in range(0,xo.__len__()):
			retval.append([xo[i],yo[i]])
		return retval

        
# TAKEN FROM MATPLOTLIB
def linspace(xmin, xmax, N):
      """
      N-element linear space between xmin and xmax

      @type xmin: float
      @param xmin: Lower boundary

      @type xmax: float
      @param xmax: Upper boundary

      @type N: int
      @param N: Number of elements

      @rtype: list
      @return: N-element list of evenly spaced points
      """
      if N==1: return [xmax]
      retval = [xmin]
      dx = (xmax-xmin)/(N-1)
      for i in range (1, N-1):
         retval += [xmin + dx*i]
      retval += [xmax]
      return retval

def extractPhi(phipsi):
      """
      Extract phi values of a string.

      @type phipsi: list
      @param phipsi: List of [phi, psi]
      """
      retval = []
      for ii in range(0, phipsi.__len__()):
         retval.append(phipsi[ii][0])
      return retval

def extractPsi(phipsi):
      """
      Extract psi values of a string.

      @type phipsi: list
      @param phipsi: List of [phi, psi]
      """
      retval = []
      for ii in range(0, phipsi.__len__()):
         retval.append(phipsi[ii][1])
      return retval

def setConstraint(phi, psi, kappa, forcefield):
    """
    Find the harmonic dihedral forces, and set their reference dihedrals

    @type phi: float
    @param phi: Reference phi value in radians

    @type psi: float
    @param psi: Reference psi value in radians

    @type args: list
    @param args: MDL force fields
    """
    flag = False
    for i in range(0, forcefield.forcetypes.__len__()):
        if (forcefield.forcetypes[i] == 'h'):
            if (not flag):
               forcefield.forcearray[i].setReferenceDihedral(phi)
	       forcefield.forcearray[i].setKappa(kappa)
               flag = True
            else:
               forcefield.forcearray[i].setReferenceDihedral(psi)
	       forcefield.forcearray[i].setKappa(kappa)


def switchPhiPsi(S):
    """
    Switch phi and psi in the passed string.
    
    @type S: list
    @param S: List of [phi, psi]
    """
    for i in range (0, len(S)):
        tmp = S[i][0]
        S[i][0] = S[i][1]
        S[i][1] = tmp

def norm(v3d):
    """
    Compute the norm of a 3-element vector.

    @type v3d: numpy.ndarray
    @param v3d: 3-element vector

    @rtype: float
    @return: The norm
    """
    return numpy.sqrt(v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2])

def norm2(v3d):
    """
    Compute the norm squared of a 3-element vector.

    @type v3d: numpy.ndarray
    @param v3d: 3-element vector

    @rtype: float
    @return: The norm squared
    """
    return v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2]

def circshift(v):
    """
    Shift all elements of the passed vector by one.
    Move the last element to the front.

    @type v: list
    @param v: Vector to shift.

    @rtype: list
    @return: Shifted vector
    """
    retval = [v[len(v)-1]]
    for i in range(0, len(v)-1):
      retval.append(v[i])
    return retval

def cumsum(v):
    """
    Return a vector with each element as the sum of all previous elements.
    Thus if [1 4 5] was passed, [1 5 10] would be returned.

    @type v: list
    @param v: Vector of elements to sum

    @rtype: list
    @return: List of cumulative sums
    """
    retval = [v[0]]
    for i in range(1, len(v)):
       retval.append(v[i]+retval[len(retval)-1])
    return retval



def M(phys, alpha, beta):
    """
    @type phys: Physical
    @param phys: The physical system

    @type alpha: int
    @param alpha: Index of the first dihedral

    @type beta: int
    @param beta: Index of the second dihedral

    @rtype: float
    @return: Value for M
    """
    #    Vector3D m = rij.cross(rkj);
    #Vector3D n = rkj.cross(rkl);
    #Vector3D fi = m * (-dVdPhi * rkj_norm / m_normsq);
    #    Vector3D fl = n * (dVdPhi * rkj_norm / n_normsq);
    #    Vector3D fj = fi * (-1 + rij_dot_rkj/rkj_normsq) - fl * (rkl_dot_rkj/rkj_normsq);
    #    Vector3D fk = - (fi + fj + fl);
    
    if (alpha == beta):
       # alpha and beta are the same dihedral, so we have four
       # terms in the sum for M
       atomI = phys.dihedral(alpha).atom1-1
       atomJ = phys.dihedral(alpha).atom2-1
       atomK = phys.dihedral(alpha).atom3-1
       atomL = phys.dihedral(alpha).atom4-1
       rij = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomI*3:atomI*3+3]
       rkj = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomK*3:atomK*3+3]
       rkl = phys.positions[atomL*3:atomL*3+3] - phys.positions[atomK*3:atomK*3+3]
       dphidxI = numpy.cross(rij,rkj)*norm(rkj)/norm2(numpy.cross(rij,rkj))
       dphidxL = -numpy.cross(rkj,rkl)*norm(rkj)/norm2(numpy.cross(rkj,rkl))
       dphidxJ = dphidxI*(-1+numpy.vdot(rij,rkj)/norm2(rkj)) - dphidxL*(numpy.vdot(rkl,rkj)/norm2(rkj))
       dphidxK = -(dphidxI + dphidxJ + dphidxL)
       
       retval = (1.0 / phys.mass(atomI+1)) * numpy.vdot(dphidxI,dphidxI)
       retval += (1.0 / phys.mass(atomJ+1)) * numpy.vdot(dphidxJ,dphidxJ)
       retval += (1.0 / phys.mass(atomK+1)) * numpy.vdot(dphidxK,dphidxK)
       retval += (1.0 / phys.mass(atomL+1)) * numpy.vdot(dphidxL,dphidxL)
       return retval
    else:
       # alpha and beta are in neighboring dihedras, so we have three
       # atoms overlapping and three terms in the sum for M
       if (alpha > beta):  # SWAP
	       tmp = alpha
	       alpha = beta
	       beta = tmp
       atomI = phys.dihedral(alpha).atom1-1
       atomJ = phys.dihedral(alpha).atom2-1
       atomK = phys.dihedral(alpha).atom3-1
       atomL = phys.dihedral(alpha).atom4-1
       rij = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomI*3:atomI*3+3]
       rkj = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomK*3:atomK*3+3]
       rkl = phys.positions[atomL*3:atomL*3+3] - phys.positions[atomK*3:atomK*3+3]
       dphidxI = numpy.cross(rij,rkj)*norm(rkj)/norm2(numpy.cross(rij,rkj))
       dphidxL = -numpy.cross(rkj,rkl)*norm(rkj)/norm2(numpy.cross(rkj,rkl))
       dphidxJ = dphidxI*(-1+numpy.vdot(rij,rkj)/norm2(rkj)) - dphidxL*(numpy.vdot(rkl,rkj)/norm2(rkj))
       dphidxK = -(dphidxI + dphidxJ + dphidxL)
       
       # DERIVATIVES FOR PSI
       atomI = phys.dihedral(beta).atom1-1
       atomJ = phys.dihedral(beta).atom2-1
       atomK = phys.dihedral(beta).atom3-1
       atomL = phys.dihedral(beta).atom4-1
       rij = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomI*3:atomI*3+3]
       rkj = phys.positions[atomJ*3:atomJ*3+3] - phys.positions[atomK*3:atomK*3+3]
       rkl = phys.positions[atomL*3:atomL*3+3] - phys.positions[atomK*3:atomK*3+3]    
       
       dpsidxI = numpy.cross(rij,rkj)*norm(rkj)/norm2(numpy.cross(rij,rkj))
       dpsidxL = -numpy.cross(rkj,rkl)*norm(rkj)/norm2(numpy.cross(rkj,rkl))
       dpsidxJ = dpsidxI*(-1+numpy.vdot(rij,rkj)/norm2(rkj)) - dpsidxL*(numpy.vdot(rkl,rkj)/norm2(rkj))
       dpsidxK = - (dpsidxI + dpsidxJ + dpsidxL)
       
       # FOR THE SUM, ONLY THE FIRST, SECOND, and THIRD ATOMS OF PSI WILL COUNT
       # THESE ARE ALSO THE SECOND, THIRD, AND FOURTH ATOMS OF PHI
       retval = 0
       retval += (1.0 / phys.mass(atomI+1)) * numpy.vdot(dphidxJ,dpsidxI)
       retval += (1.0 / phys.mass(atomJ+1)) * numpy.vdot(dphidxK,dpsidxJ)
       retval += (1.0 / phys.mass(atomK+1)) * numpy.vdot(dphidxL,dpsidxK)
       return retval	
    

# TREVOR'S SMOOTHING FUNCTION
def L(z, p):
   result = 0
   if (p == 1):
      return result
   for q in range(2, p+1):
      result += numpy.sqrt((z[q-1][0]-z[q-1-1][0])**2+(z[q-1][1]-z[q-1-1][1])**2)
   return result

def l(z, p):
   R = len(z)
   return (p-1)*L(z, R)/(R-1)

def q(z, p):
   tmpQ = 2
   myl = l(z, p)
   while (not (L(z, tmpQ-1) < myl and myl <= L(z, tmpQ))):
         tmpQ += 1
   return tmpQ

def reparamTrevor(z):
   S = []
   S.append(z[0])
   for p in range(2, len(z)):
      myQ = q(z, p)
      newphi = z[myQ-1-1][0] + (l(z, p)-L(z, myQ-1))*(z[myQ-1][0]-z[myQ-1-1][0])/numpy.sqrt((z[myQ-1][0]-z[myQ-1-1][0])**2+(z[myQ-1][1]-z[myQ-1-1][1])**2)
      newpsi = z[myQ-1-1][1] + (l(z, p)-L(z, myQ-1))*(z[myQ-1][1]-z[myQ-1-1][1])/numpy.sqrt((z[myQ-1][0]-z[myQ-1-1][0])**2+(z[myQ-1][1]-z[myQ-1-1][1])**2)
      S.append([newphi, newpsi])
   S.append(z[len(z)-1])
   return S

# CHRIS' SMOOTHING FUNCTION
def reparam(z):
    """
    High level smoothing/reparameterization routine.  Invokes spline functions.

    @type x: list
    @param x: X values

    @type y: list
    @param y: Y values

    @rtype: list
    @return: Smoothed and reparameterized list of [x, y] pairs
    """
    x = extractPhi(z)
    y = extractPsi(z)
    dx = []
    dy = []
    cx = circshift(x)
    cy = circshift(y)
    for i in range(0, len(x)):
       dx.append(x[i]-cx[i])
       dy.append(y[i]-cy[i])
    dx[0] = 0
    dy[0] = 0

    dist_vec = []
    for i in range(0, len(dx)):
       dist_vec.append(math.sqrt(dx[i]*dx[i]+dy[i]*dy[i]))

    lxy = cumsum(dist_vec)
    for i in range(0, len(lxy)):
       lxy[i] /= lxy[len(lxy)-1]

    mySpline = SplineInterpolator()
    S_x = mySpline.Spline(lxy, x)
    g1 = linspace(lxy[0], lxy[len(lxy)-1], 100)
    vx = mySpline.EvalSpline(lxy, x, S_x, g1)
    S_x = mySpline.SmoothingAndReparameterization(lxy, x, g1, vx, len(lxy))

    S_y = mySpline.Spline(lxy, y)
    vy = mySpline.EvalSpline(lxy, y, S_y, g1)
    S_y = mySpline.SmoothingAndReparameterization(lxy, y, g1, vy, len(lxy))

    S = []
    for i in range(0, len(S_x)):
       S.append([S_x[i][1], S_y[i][1]])
    return S
    
