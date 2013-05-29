import warnings
warnings.filterwarnings(action='ignore',
                        message='.*has API version.*',
                        category=RuntimeWarning)

import Forces
import Constants
import Vector3DBlock
import PARReader
import PSFReader
import PDBReader
import XYZReader
import EigenvectorReader
import MathUtilities
import TopologyUtilities
import GenericTopology
#import _TopologyUtilities
import _GenericTopology
import numpy
import sys
import ProtoMolApp
def deepcopy(a, b):
   if (len(b) != len(a)): # if the length of a is not equal to the length of a
      b.resize(len(a)) # resize the length of a to the size of b
   for i in range (0, len(a)): # for-loop i from 0 to the length of a
      b[i] = a[i] # set a[i] equal to b[i]


class MDVec(numpy.ndarray): # takes a paramater of an ndarray
   def __eq__(self, b):
    if (len(self) != len(b)): # if the length of self is not equivalent to the length of b
      self.resize(len(b)) # resize the length of b to the resize self
    for i in range (0, len(b)): # for-loop from 0 to the length of b
      self[i] = b[i] # set b[i] equivalent to self[i]
      

class Atom:
   """
   An atom in the system.
   """
   def __init__(self, n, s, rs, rn, an, at, c, m):
      self.number = n   #: Atom number
      self.seg_id = s   #: Segment identifier
      self.residue_sequence = rs  #: Residue sequence
      self.residue_name = rn     #: Residue name
      self.atom_name = an        #: Atom name
      self.atom_type = at        #: Atom type
      self.charge = c            #: Charge [e]
      self.mass = m              #: Mass [amu]

class Bond:
   """
   A two-atom bond.
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2

class Angle:
   """
   A three-atom angle.
   """
   def __init__(self, n, a1, a2, a3):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3

class Dihedral:
   """
   A four-atom dihedral.
   """
   def __init__(self, n, a1, a2, a3, a4):
      self.number = n    #: Dihedral number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3
      self.atom4 = a4    #: Atom index 4

class Improper:
   """
   A four-atom improper.
   """
   def __init__(self, n, a1, a2, a3, a4):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3
      self.atom4 = a4    #: Atom index 4

class HDonor:
   """
   An H+ donor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Donor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2

class HAcceptor:
   """
   An H+ acceptor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Acceptor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2


global apps
apps = []     

class Physical:
   """
   Defines a physical system: positions, velocities, temperature,
   boundary conditions, etc.
   """


   def __init__(self):
      #####################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.seed = 1234               #: Random number seed
      self.exclude = "1-4"           #: Exclusion Pairs
      self.cellsize = 6.5            #: Cell size
      self.bc = "Periodic"           #: Boundary conditions
      self.remcom = 0                #: Remove COM motion?
      self.remang = 0                #: Remove angular momentum?
      self.defaultCBV = True         #: Use default cell basis vectors (PBC)
      self.time = 0                  #: Current time
      self.cB1 = numpy.ndarray(3)    #: Cell basis vector 1
      self.cB1.fill(0)
      self.cB2 = numpy.ndarray(3)    #: Cell basis vector 2
      self.cB2.fill(0)
      self.cB3 = numpy.ndarray(3)    #: Cell basis vector 3
      self.cB3.fill(0)
      self.cO = numpy.ndarray(3)     #: Cell origin
      self.cO.fill(0)
      self.temperature = 300         #: Kelvin temperature
      self.masses = numpy.ndarray(0)      #: Diagonal mass matrix
      self.invmasses = numpy.ndarray(0)   #: Diagonal inverse mass matrix
      self.masssum = 0                    #: Sum over all atomic masses
      self.accept = 0
      self.reject = 0
      ProtoMolApp.ProtoMolApp.turnOffHints()
      # NOTE: self.positions and self.velocities are also
      # available numpy arrays.
      #####################################################################

      ############################################
      # NOTE: USER SHOULD NOT TOUCH THESE!
      # THESE ARE NUMPY ARRAY WRAPPERS
      # USER SHOULD ACCESS THROUGH
      # self.positions and self.velocities
      self.__dict__['myTop'] = GenericTopology.T_Periodic() # 'mytop' is a GenericTopology.T_Periodic()
      self.__dict__['posvec'] = Vector3DBlock.Vector3DBlock() # 'posvec' is  a Vector3DBlock.Vector3DBlock()
      self.__dict__['velvec'] = Vector3DBlock.Vector3DBlock() # 'velvec' is a Vector3DBlock.Vector3DBlock()
      self.__dict__['myPAR'] = PARReader.PAR() # 'myPAR' is a PARReader.PAR()
      self.__dict__['myPSF'] = PSFReader.PSF() # 'myPSF' is a PSFReader.PSF()
      self.__dict__['myPDB'] = PDBReader.PDB() # 'myPDB' is a PDBReader.PDB()
      self.__dict__['myXYZ'] = XYZReader.XYZ() # 'myXYZ' is a XYZReader.XYZ()
      self.__dict__['myEig'] = EigenvectorReader.EigenvectorInfo() # myEig' is a EigenvectorReader.EigenvectorInfo()
      ############################################

      #self.positions = numpy.ndarray(0)   #: Atomic position vector
      #self.velocities = numpy.ndarray(0)  #: Atomic velocity vector

      self.dirty = 1   #: Dirty bit
   
   def __del__(self):
      pass 

   def __str__(self):
      return "Physical.Physical"

   def __repr__(self):
      return "Physical.Physical"

   # Copy which avoids object assignment
   def copy(self, retval, forces="", dt=-1):
      """
      Perform a deep copy, avoid reference assignment
      """
      #retval = Physical()
      retval.seed = self.seed # set seed to retval.seed
      retval.exclude = self.exclude # set exclude to retval.exclude
      retval.cellsize = self.cellsize # set cellsize to retval.cellsize
      retval.bc = self.bc # set bc to retval.bc
      retval.defaultCBV = self.defaultCBV # set defaultCBV to retval.defaultCBV
      retval.remcom = -1 # set -1 to retval.remcom
      retval.remang = -1 # set -1 to retval.remang
      retval.defaultCBV = self.defaultCBV #set defaultCBV to retval.defaultCBV
      retval.myTop = self.myTop # set myTop to retval.myTop
      #if (retval.bc == "Periodic"):
      #   retval.myTop = GenericTopology.T_Periodic()
      #else:
      #   retval.myTop = GenericTopology.T_Vacuum()
      retval.cB1 = self.cB1.copy() # copy cB1 and set it to retval.cB1
      retval.cB2 = self.cB2.copy() # copy cB2 and set it to retval.cB2
      retval.cB3 = self.cB3.copy() # copy cB3 and set it to retval.cB3
      retval.cO = self.cO.copy()   # copy cO and set it to retval.cO
      retval.temperature = self.temperature # set temperature to retval.temperature
      retval.masses = self.masses.copy() # copy masses and set it to retval.masses
      retval.invmasses = self.invmasses.copy() # copy invmasses and set it to retval.invmasses
      retval.masssum = self.masssum # set masssum to retval.masssum
      retval.posvec.resize(self.posvec.size()) # resize 
      for i in range(len(self.positions)): # for-loop i to the length of self.positions
         retval.positions[i] = self.positions[i] # set self.positions[i] equal to retval.positions[i]
      retval.velvec.resize(self.velvec.size()) # resize
      for i in range(len(self.velocities)): # for-loop i to the length of self.velocities
         retval.velocities[i] = self.velocities[i] # set self.velococities[i] equal to retval.velocities[i]
      retval.myPAR = self.myPAR # set myPAR to retval.myPAR
      retval.myPSF = self.myPSF # set myPSF to retval.myPSF
      retval.myPDB = self.myPDB # set myPDB to retval.myPDB
      retval.myEig = self.myEig # set myEig to retval.myEig
      retval.accept = self.accept # set accept to retval.accept
      retval.reject = self.reject # set reject to retval.reject

      if (dt != -1): # if dt does not equal -1
         retval.app = ProtoMolApp.ProtoMolApp() # set ProtoMolApp to retval.app
         f = Forces.Forces() # forces are set to f
      	 retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, f.energies, dt) # makeApp
         apps.append(retval.app) # append retval.app to apps
      #retval.build()
      #return retval
      

   def deepcopy(self, forces, dt):
         retval = self.copy() # set the copy of self to retval
         retval.app = ProtoMolApp.ProtoMolApp() # set ProtoMolApp to retval.app
         retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, forces.energies, dt) # makeApp

      
   # SPECIAL ACCESSOR FOR self.positions or self.velocities
   # TO GET DATA FROM WRAPPERS
   def __getattr__(self, name):
      if (name == 'positions'): # if the name = 'positions'
	 return (self.__dict__['posvec'].getC()).view(MDVec) # return the position vector
      elif (name == 'velocities'):# if the name = 'velocities'
         return (self.__dict__['velvec'].getC()).view(MDVec)  # return the velocity vector
      elif (name == 'time'): # if the name = 'time'
         return self.myTop.time # return myTop time
      else: # if name is not set to 'positions', 'velocities'  or 'time' then 
         print name # print name
         return self.__dict__[name] # return  __dict_[name]

   # SPECIAL ASSIGNMENT FOR self.positions or self.velocities
   # TO SET DATA IN WRAPPERS
   def __setattr__(self, name, val):
      firsttime = False # firsttime is set to False
      if (not self.__dict__.has_key(name)): # if self.__dict__ does not have key(name) 
         firsttime = True # firsttime is set to True
      if (name == 'positions'): # if name = 'positions'
	 for i in range (0, len(self.positions)): # loop from 0 to the length of positions
	    self.positions[i] = val[i] # set val[i] to positions[i]
	 #self.__dict__['posvec'].setC(val)
      elif (name == 'velocities'): # if name = 'velocities'
	 for i in range (0, len(self.velocities)): # loop from 0 to length of velocities
	    self.velocities[i] = val[i] # set val[i] to velocities[i]
	 #self.__dict__['velvec'].setC(val)
      elif (name == 'bc'): # if the name = 'bc'
         self.__dict__['bc'] = val         # then set val to self.__dict__['bc']
         if (val == "Vacuum"): # is val = "Vacuum" 
           self.myTop = GenericTopology.T_Vacuum() # then set GenericTopology.T_Vacuum() to myTop
         else: # if val is not a "Vacuum"
           self.myTop = GenericTopology.T_Periodic() # set GenericTopology.T_Periodic() to myTop
         if (not firsttime): # is this is not firsttime
            self.build() # build()
      elif (name == 'time'): # if name = 'time'
         val /= Constants.invTimeFactor() # val = val/~.0205
         self.myTop.time = val # set val to myTop time
      elif (name == 'seed' or  # is name is 'seed','exclude','cellsize','remcom', or 'remang'
            name == 'exclude' or
            name == 'cellsize' or
            name == 'remcom' or
            name == 'remang'):
         self.__dict__[name] = val # set val equal to the name
         if (not firsttime): # if not firsttime
            self.build() # build
      else: # any other name
         self.__dict__[name] = val # set value to name

   # RESETS SIMULATION STATE TO DEFAULTS
   # ALSO RESETS TIME TO ZERO
   def reset(self):
      """
      Reset all member variables to default values
      """
      self.__dict__['seed'] = 1234               # Random number seed
      self.__dict__['exclude'] = "1-4"           # Exclusion Pairs
      self.__dict__['cellsize'] = 6.5            # Cell size
      self.__dict__['bc'] = "Periodic"           # Boundary conditions
      self.__dict__['remcom'] = 0                # Remove COM motion?
      self.__dict__['remang'] = 0                # Remove angular momentum?
      self.__dict__['defaultCBV'] = True         # Use default cell basis vectors (PBC)
      self.__dict__['myTop'] = Topology.T_Periodic()
      self.time = 0 # set time to 0
      self.cB1.fill(0) # fill cB1 with 0
      self.cB2.fill(0) # fill cB2 with 0
      self.cB3.fill(0) # fill cB3 with 0
      self.cO.fill(0)  # fill cO with 0
      self.__dict__['temperature'] = 300         # Kelvin temperature

   # SYSTEM PRESSURE.
   def pressure(self, forces):
      """
      Pressure of the system.

      @type forces: Forces
      @param forces: MDL Forces object

      @rtype: float
      @return: System pressure
      """
      TopologyUtilities.addVelocityVirial(forces.energies, self.myTop, self.velvec) # addVelocity Virial to TopologyUtilities
      return TopologyUtilities.velocityVirial(self.myTop, self.velvec).pressure(self.myTop.getVolume(self.posvec)) # return system pressure

   # SYSTEM VOLUME (AA^3)
   def volume(self):
      """
      Volume of the system.

      @rtype: float
      @return: System volume
      """
      return self.myTop.getVolume(self.posvec) # return system volume


   # SIZE OF THE SYSTEM
   def N(self):
      """
      Number of atoms.

      @rtype: int
      @return: Number of atoms.
      """
      return self.numAtoms() # return Number of atoms
   
   def numAtoms(self):
      """
      Number of atoms.

      @rtype: int
      @return: Number of atoms.
      """
      return self.myPSF.numAtoms() # return Number of atoms

   # NUMBER OF TWO-ATOM BONDS.
   def numBonds(self):
      """
      Number of bonds.

      @rtype: int
      @return: Number of bonds
      """
      return self.myPSF.numBonds() # return Number of bonds

   # NUMBER OF THREE-ATOM ANGLES.
   def numAngles(self):
      """
      Number of angles

      @rtype: int
      @return: Number of angles
      """
      return self.myPSF.numAngles() # return number of angles

   # NUMBER OF FOUR-ATOM DIHEDRALS.
   def numDihedrals(self):
      """
      Number of dihedrals

      @rtype: int
      @return: Number of dihedrals
      """      
      return self.myPSF.numDihedrals() # return Number of dihedrals

   # NUMBER OF FOUR-ATOM IMPROPERS.
   def numImpropers(self):
      """
      Number of impropers

      @rtype: int
      @return: Number of impropers
      """
      return self.myPSF.numImpropers() # return Number of impropers

   # NUMBER OF H-BOND DONORS.
   def numDonors(self):
      """
      Number of hydrogen donors (for H+ bonding)

      @rtype: int
      @return: Number of hydrogen donors
      """
      return self.myPSF.numDonors()# return Number of hydrogen donors

   # NUMBER OF H-BOND ACCEPTORS.
   def numAcceptors(self):
      """
      Number of hydrogen acceptors (for H+ bonding)

      @rtype: int
      @return: Number of hydrogen acceptors
      """
      return self.myPSF.numAcceptors() # return Number of hydrogen acceptors

   # GET ATOM #(index)
   def atom(self, index):
      """
      Get an atom at the passed index.

      @type index: int
      @param index: Atom index (1 to N)

      @rtype: Atom
      @return: The atom at the passed index.
      """
      return Atom(self.myPSF.getAttributeInt("atom", index-1, "number"), # get atom at the passed index, which consits of number
                  self.myPSF.getAttributeString("atom", index-1, "seg_id"), #segmented id
                  self.myPSF.getAttributeInt("atom", index-1, "residue_sequence"), #residue sequence
                  self.myPSF.getAttributeString("atom", index-1, "residue_name"),# residue name
                  self.myPSF.getAttributeString("atom", index-1, "atom_name"),# atom name
                  self.myPSF.getAttributeString("atom", index-1, "atom_type"),# atom type
                  self.myPSF.getAttributeReal("atom", index-1, "charge"),# charge
                  self.myPSF.getAttributeReal("atom", index-1, "mass"))# and mass


   # GET BOND #(index)
   def getBond(self, index):
      """
      Get a bond at the passed index.

      @type index: int
      @param index: Bond index

      @rtype: Bond
      @return: The bond at the passed index.
      """
      return Bond(self.myPSF.getAttributeInt("bond", index-1, "number"), # get bond at passed index including number
                  self.myPSF.getAttributeInt("bond", index-1, "atom1"),  # atom1
                  self.myPSF.getAttributeInt("bond", index-1, "atom2"))  # and atom2

   # GET ANGLE #(index)
   def getAngle(self, index):
      """
      Get an angle at the passed index.

      @type index: int
      @param index: Angle index

      @rtype: Angle
      @return: The angle at the passed index.
      """
      return Angle(self.myPSF.getAttributeInt("angle", index-1, "number"), # get the angle at the passed index including number
                   self.myPSF.getAttributeInt("angle", index-1, "atom1"),  # atom1
                   self.myPSF.getAttributeInt("angle", index-1, "atom2"),  # atom2
                   self.myPSF.getAttributeInt("angle", index-1, "atom3"))  # and atom3

   # GET DIHEDRAL #(index)
   def getDihedral(self, index):
      """
      Get a dihedral at the passed index.

      @type index: int
      @param index: Dihedral index

      @rtype: Dihedral
      @return: The dihedral at the passed index.
      """
      return Dihedral(self.myPSF.getAttributeInt("dihedral", index-1, "number"),# get the dihedral at the passed index including number
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom1"), # atom1
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom2"), # atom2
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom3"), # atom3
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom4")) # atom4

   # GET IMPROPER #(index)
   def getImproper(self, index):
      """
      Get an improper at the passed index.

      @type index: int
      @param index: Improper index

      @rtype: Improper
      @return: The improper at the passed index.
      """
      return Improper(self.myPSF.getAttributeInt("improper", index-1, "number"),# get the improper at the passed index including number
                      self.myPSF.getAttributeInt("improper", index-1, "atom1"), # atom1
                      self.myPSF.getAttributeInt("improper", index-1, "atom2"), # atom2
                      self.myPSF.getAttributeInt("improper", index-1, "atom3"), # atom3
                      self.myPSF.getAttributeInt("improper", index-1, "atom4")) # atom4


   # GET DONOR #(index)
   def getDonor(self, index):
      """
      Get an H+ donor at the passed index.

      @type index: int
      @param index: H+ donor index

      @rtype: HDonor
      @return: The H+ donor at the passed index.
      """
      return HDonor(self.myPSF.getAttributeInt("donor", index-1, "number"),# get the donor at the passed index including number
                    self.myPSF.getAttributeInt("donor", index-1, "atom1"), # atom1
                    self.myPSF.getAttributeInt("donor", index-1, "atom2")) # and atom2

   # GET ACCEPTOR #(index)
   def getAcceptor(self, index):
      """
      Get an H+ acceptor at the passed index.

      @type index: int
      @param index: H+ acceptor index

      @rtype: HAcceptor
      @return: The H+ acceptor at the passed index.
      """
      return HAcceptor(self.myPSF.getAttributeInt("acceptor", index-1, "number"),# get the acceptor at the passed index including number
                       self.myPSF.getAttributeInt("acceptor", index-1, "atom1"), # atom 1
                       self.myPSF.getAttributeInt("acceptor", index-1, "atom2")) # and atom2



   # GET THE MASS (AMU) OF A SPECIFIC ATOM
   def mass(self, atom):
      """
      Mass of an atom.

      @type atom: int
      @param atom: Atom index (1 to N)

      @rtype: float
      @return: Atom mass [amu]
      """
      return self.myPSF.getAttributeReal("atom", atom-1, "mass") # return the Atomic mass

   def charge(self, atom):
      """
      Charge of an atom.

      @type atom: int
      @param atom: Atom index (1 to N)

      @rtype: float
      @return: Atomic charge [e]
      """
      return self.myPSF.getAttributeReal("atom", atom-1, "charge") # return the Atomic charge [e]

   # SYSTEM TEMPERATURE.
   def temperature(self):
      """
      System temperature (K)

      @rtype: float
      @return: Kelvin temperature
      """
      return TopologyUtilities.temperature(self.myTop, self.velvec) # return the Kelvin temperature

   def dihedral(self, index):
      """
      Dihedral angle (rad, -PI to PI) at passed index

      @type index: int
      @param index: Dihedral index
      
      @rtype: float
      @return: Dihedral angle in radians
      """
      myPhi = TopologyUtilities.computePhiDihedral(self.myTop, self.posvec, index-1) # computePhiDihedra of ToplogyUtilities and set it to myPhi
      if (myPhi > numpy.pi): # if myPhi is greater than numpy.pi
         myPhi -= 2*numpy.pi # myPhi = myPhi - 2*numpy.pi
      elif (myPhi < -numpy.pi): # if myPhi is less than the -numpy.pi
         myPhi += 2*numpy.pi # myPhi = myPhi + 2*numpy.pi
      return myPhi # return myPhi


   def randomVelocity(self, T):
      """
      Assign random velocities.

      @type T: float
      @param T: Kelvin temperature.
      """
      TopologyUtilities.randomVelocity(T, self.myTop, self.velvec, self.seed) # get the randomVelocity of TopologyUtilities

      
   def updateCOM_Momenta(self):
      """
      Update center of mass and angular momentum
      """
      TopologyUtilities.buildMolecularCenterOfMass(self.posvec,self.myTop) # update center of mass
      TopologyUtilities.buildMolecularMomentum(self.velvec,self.myTop) #update angular momentum

   def build(self):
      """
      Build the physical data.
      """
      tm = -1 # tm is set to -1
      if (hasattr(self, "myTop")): # if self has the attribute "myTop"
          tm = self.myTop.time # tm = myTop.time
      if (self.bc == "Periodic"): # if bc is equal to "Periodic"
          self.myTop.setCellSize(self.cellsize) # setCellSize to cellsize
          if (self.defaultCBV): # if self is set to defacultCBV
             # TEMPORARY STRUCTURES USED TO COMPUTE
             # BOUNDING BOX
             v1 = numpy.ndarray(3) # v1,v2,v3,and v4 are numpy.ndarray of size 3
             v2 = numpy.ndarray(3)
             v3 = numpy.ndarray(3)
             v4 = numpy.ndarray(3)
             v1.fill(sys.maxint) # fill v1 with maxint
             v2.fill(-sys.maxint)# fill v2 with the negative of the maxint
             # BOUNDING BOX
             i = 0 # i set to 0
             while (i < numpy.size(self.positions)): #while i is less than the size of positions
                if (self.positions[i] < v1[0]): # and if self.positions[i] is less than v1[0]
                   v1[0] = float(self.positions[i]) # then v1[0] = self.positions[i]
                if (self.positions[i] > v2[0]): # and if self.positions[i] greater than v2[0]
                   v2[0] = float(self.positions[i]) # then v2[0] = self.positions[i]
                if (self.positions[i+1] < v1[1]): # and if self.positions[i+1] is less than v1[1] 
                   v1[1] = float(self.positions[i+1]) # then v1[1] = positions[i+1]
                if (self.positions[i+1] > v2[1]): # and if self.positions[i+1] is greater than v2[1]
                   v2[1] = float(self.positions[i+1]) # then v2[1] = positions[i+1]
                if (self.positions[i+2] < v1[2]): # and if self.positions[i+2] is less than v1[1]
                   v1[2] = float(self.positions[i+2]) # then v1[2] = positions[i+2]
                if (self.positions[i+2] > v2[2]): # and if self.positions[i+2] is greater than v2[1]
                   v2[2] = float(self.positions[i+2]) # then v2[2] = positions[i+2]
                i += 3 # iterate i by 3
             v4.fill(Constants.periodicBoundaryTolerance()/2.0) # fill v4 with periodicBoundaryTolerance() /2.0
             v1 = v1 - v4 # set v1 - v4 to v1
             v2 = v2 + v4 # set v2 + v4 to v2
             v3 = v2 - v1 # set v2 - v1 to v3
             self.cB1[0] = v3[0] # set v3[0] to cB1[0] 
             self.cB2[1] = v3[1] # set v3[1] to cB2[1]
             self.cB3[2] = v3[2] # set v3[2] to cB3[2]
             self.cO = v1 + v3 * 0.5 # set v1 + v3 *.5 to cO
             self.myTop.setBC(self.cB1[0],self.cB1[1],self.cB1[2],self.cB2[0],self.cB2[1],self.cB2[2],self.cB3[0],self.cB3[1],self.cB3[2],self.cO[0],self.cO[1],self.cO[2]) # setBC
      self.myTop.setExclusion(self.exclude) # setExclusion
      if (self.myPSF.numAtoms() > 0 and hasattr(self.myPAR, 'readFlag')): # if number of Atoms is greater than 0 and myPAR has 'redFlag' attribute
	 GenericTopology.buildTopology(self.myTop, self.myPSF, self.myPAR, 0, self.myTop.makeSCPISM()) # buildTopology
	 if (numpy.size(self.velocities) == 0):  if the size of self.velocities = 0 # if the numpy size of velocities is equal to 0 and 
            # if the size of velocities = 0

            # This function actually returns a value, since it is
            # unused the runtime environment detects a memory leak.
            # disown() removes this concern.
            MathUtilities.randomNumber(self.seed) 
            # NOTE TO SELF: THIS WAS COMMENTED OUT FOR PM 3
            # IT MAY EFFECT RANDOM NUMBER CONSISTENCY
            #aaa = MathUtilities.randomNumberFirst(self.seed, 1)
	    TopologyUtilities.randomVelocity(self.temperature, self.myTop, self.velvec, self.seed) # TopologyUtilities randomVelocity
	 if (self.remcom >= 0): # if remcom is greater than or equal to 0
	   TopologyUtilities.removeLinearMomentum(self.velvec, self.myTop).disown() # removeLinearMomentum
         if (self.remang >= 0): # if remang is greater than or equal to 0
	   TopologyUtilities.removeAngularMomentum(self.posvec, self.velvec, self.myTop).disown() # removeAngularMomentum
      #if (self.bc == "Periodic"):
      #     self.myTop.setCellSize(self.cellsize)
      # COMPUTE INV MASS MATRIX
      temp = list() # list() is set to temp
      ii = 0 
      while ii < self.numAtoms()*3: # while-loop from ii less than number of Atoms *3 
          temp.append(0.0) # append 0.0 to temp
          ii += 1 # iterate ii by 1
      ii = 0
      self.invmasses.resize(self.numAtoms()*3) # resize invmasses to the number of Atoms multiplied by 3N

      # iteration over number of atoms *3 = 3N
      # we get the atom by dividing by 3 since there are 3 parts to the atoms {x,y,z}
      # F = Ma
      # where inverse mass is 3N by 3N
      # force is 3N by 1 and so is acceleration
      while ii < self.numAtoms()*3: # going from 0 to 3N
          im = 1.0 / self.myPSF.getMass(ii/3) 
          self.invmasses[ii] = im     # invmasses 3n x 3n
          self.invmasses[ii+1] = im
          self.invmasses[ii+2] = im
          ii += 3 # iterating by 3

      # COMPUTE MASS MATRIX
      temp = list() # set list equivalent to temp
      ii = 0
      while ii < self.numAtoms()*3:
          temp.append(0.0) # append 0.0 to temp
          ii += 1 # iterate ii by 1
      ii = 0
      self.masses.resize(self.numAtoms()*3)
      while ii < self.numAtoms()*3: # compare the iteration with the number of Atoms *3
          m = self.myPSF.getMass(ii/3) # m is equal to the myPSF mass over a period of iteration
          temp[ii] = m # set m to temp[ii]
          self.masses[ii] = temp[ii] # set temp[ii] to masses[ii]
          temp[ii] = 0.0 # set temp[ii] to 0
          temp[ii+1] = m # set m to temp[ii+1]
          self.masses[ii+1] = temp[ii+1] # set temp[ii+1] to masses[ii+1]
          temp[ii+1] = 0.0 # set temp[ii+1] to 0.0
          temp[ii+2] = m # set m to temp[ii+2]
          self.masses[ii+2] = temp[ii+2] # set temp[ii+2] to masses[ii+2]
          temp[ii+2] = 0.0 # set the temp[ii+2] to 0.0
          ii += 3 # iterate ii by 3


      # COMPUTE MASS SUM
      self.masssum = 0 # set masssum to 0
      ii = 0 # set ii to 0
      while ii < self.numAtoms(): # compare iteration with number of Atoms.
         self.masssum += self.myPSF.getMass(ii) # masssum = masssum + myPSF.getMass(ii) 
         ii += 1 # iterate ii

      # SET SELFICAL TIME
      if (tm != -1): # if tm does not equal -1
          self.myTop.time = tm # set tm to myTop time

      self.dirty = 0  # clean now!


