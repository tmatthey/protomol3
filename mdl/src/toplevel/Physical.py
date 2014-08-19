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
#from GPUMatrix import GpuMatFunctions, GpuMatrix
#from theano import config, shared

# This is a deepcopy function for two Python lists a and b
# If you simply set b = a that copies the reference and if you change a you change b
# In this situation, we do an element-by-element copy from a to b which
# is desirable sometimes.
def deepcopy(a, b):
   if (len(b) != len(a)): # if the length of a is not equal to the length of a
      b.resize(len(a)) # resize the length of a to the size of b
   for i in range (0, len(a)): # for-loop i from 0 to the length of a
      b[i] = a[i] # set a[i] equal to b[i]


# The MDVec class is almost identical to numpy.ndarray
# Note that it inherits from numpy.ndarray
# However, we add the __eq__ function to be able to control what happens
# when we set a = b and a and b are both MDVecs.
# This calls a.__eq__(b) behind the scenes
# And our function does an element-by-element copy from a to b
# Otherwise it will set the references equal and we have the problem above.
class MDVec(numpy.ndarray): # takes a paramater of an ndarray
   def __eq__(self, b):
    if (len(self) != len(b)): # if the length of self is not equivalent to the length of b
      self.resize(len(b)) # resize the length of b to the resize self
    for i in range (0, len(b)): # for-loop from 0 to the length of b
      self[i] = b[i] # set b[i] equivalent to self[i]
   def toNumpy(self,MDVec):
      self = numpy.array(MDVec) 	   
   def toGPUMAtrix(self,MDVec):
      self = numpy.array(MDVec)
      self = gpu.GPUMatrix(self)

# Information about a single atom in our system.
# 1. The atom number is a unique identifier, starting from zero and going until N-1
# where N is the number of atoms in our system.
# 2. seg_id is the segment identifier, which for systems of one molecule will be the same
# for all atoms.  But we may incorporate explicit solvent and have many water molecules,
# in this case you will have different seg_id values for the solvent and solute.
# 3. The residue_sequence and residue_number define the residue to which this atom belongs.
#    A protein molecule is made up of residues (or amino acids) of the following form:
#    H       H      O 
#    |       |      ||
#    N-------C------C--------OH
#    |       |
#    H       R
#
# R is what differentiates amino acids from each other.  There are twenty total:
# Glycine, Alanine, Valine, Leucine, Isoleucine, Methionine, Tryptophan,
# Serine, Threonine, Cysteine, Tyrosine, Asparagine, Glutamine, Aspartic Acid,
# Glutamic Acid, Lysine, Arginine, Histidine, Phenylalanine, and Proline.
# Proteins are long sequences of amino acids combined through dehydration synthesis.
# 4. Atom_name is the name of the atom, i.e. C ("Carbon")
# 5. Atom_type is the type of the atom - sometimes this is simply equal to the name
# but other times it is a specific atom, for example the first carbon above is termed
# an alpha-carbon (CA) because it is bound to the residue R.
# 6. Charge is the atomic charge.
# 7. Mass is the atomic mass.
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


# A bond in our system consists of:
# 1. A unique identifier, starting from zero
# 2. The atom numbers of the two atoms involved in the bond.
#    These correspond to the 'number' data member of the Atom class above.
class Bond:
   """
   A two-atom bond.
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2


# An angle in our system consists of:
# 1. A unique identifier, starting from zero
# 2. The atom numbers of the three atoms involved in the angle.
#    These correspond to the 'number' data member of the Atom class above.
class Angle:
   """
   A three-atom angle.
   """
   def __init__(self, n, a1, a2, a3):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3


# A dihedral in our system consists of:
# 1. A unique identifier, starting from zero
# 2. The atom numbers of the four atoms involved in the dihedral.
#    These correspond to the 'number' data member of the Atom class above.
# The structure of a dihedral is: A---B---C---D
# The dihedral angle is the angle between the plane formed by A---B---C
# and the plane formed by B---C---D
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


# An improper in our system consists of:
# 1. A unique identifier, starting from zero
# 2. The atom numbers of the four atoms involved in the improper.
#    These correspond to the 'number' data member of the Atom class above.
# The four atoms in an improper are bonded differently compared to dihedrals.
# The structure is:
#        C
#        |
#    B---A---D
#
# Thus there is a central atom bound to all three
# The improper angle is thus now the angle between the plane formed by A---B---C
# and the plane formed by A---C---D   (not B---C--D as above)
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


# Hydrogen bonding between non-sequential amino acids of a protein
# between the residue (R, see above) of a polar amino acid (e.g.
# Serine, Threonine, Cysteine, Tyrosine, Asparagine and Glutamate) and
# those that have an exposed negatively charged oxygen in their residue
# (e.g., Aspartic Acid, Glutamic Acid, Lysine, Arginine and Histidine).
# When non-sequential amino acids bind in this way they can form structures
# like alpha-helices and beta-sheets.
# In the case of an HDonor, we have:
# 1. A unique identifier starting from zero (number)
# 2. atom1, the index of the H+ ion (corresponding to the atom numbers in the Atom class)
# 3. atom2, the atom to which the H+ is bonded (the donor)
class HDonor:
   """
   An H+ donor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Donor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2


# In the case of an acceptor, we have:
# 1. A unique identifier starting from zero (number)
# 2. atom1, the index of the O- ion (corresponding to the atom numbers in the Atom class)
# 3. atom2, the atom to which the O- is bonded (the acceptor).
class HAcceptor:
   """
   An H+ acceptor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Acceptor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2


# This structure holds objects of type ProtoMolApp 
# This is the central SWIG-wrapped structure from ProtoMol, which contains
# the core of all the information used to run a simulation
# We have to make sure that these structures are not deallocated until
# the Python interpreter is closed, therefore we keep them global.
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
      
      # Random number seed (all simulations with this same seed will result
      # in identical results for output, even if they have random processes
      # internally.
      self.seed = 1234               #: Random number seed
      self.gpu = False
      # Possible values for this include: "1-2", "1-3", "1-4" and "scaled1-4"
      # If set to 1-2, we exclude (do not compute) electrostatic interactions between atoms directly 
      # bonded to each other (because we took this into account with the bond force)
      # If set to 1-3, we exclude 1-2 plus electrostatic interactions between atoms
      # in the same angle (because we have an Angle force)
      # If set to 1-4, we exclude 1-3 plus electrostatic interactions between atoms
      # in the same dihedral or improper (because we have those forces).  This is the default.
      # If set to scaled1-4, we exclude 1-3 and scale electrostatic forces between atoms
      # in the same dihedral or improper by a small value (but do not completely neglect them).
      self.exclude = "1-4"           #: Exclusion Pairs

      ############################################################################################
      # The following parameters are only used for Periodic Boundary Conditions (PBC)
      # PBC is mainly used within solvated systems so that we can achieve wraparound in 3D
      # If we envision a solvated system in a 3D 'box' or 'cell' for instance of size 6.5 Angstroms
      # (the default value for the cell size) -- a water molecule 0.5 Angstroms from the 'top'
      # of this box (6) should interact with a water molecule 1 Angstrom from the 'bottom' of the box
      # assuming a distance of 1.5 Angstroms between them (not 5 Angstroms) because of wraparound 
      # in a closed system.
      self.cellsize = 6.5            #: Cell size
      self.bc = "Periodic"           #: Boundary conditions
      
      # Cell basis vector information
      # These define the structure of the cells (default is cubic)
      # Initially we fill these three vectors as (0,0,0) but they will be set to their 
      # default values in build().  Otherwise you MUST set them yourself.
      self.defaultCBV = True         #: Use default cell basis vectors (PBC)
      self.cB1 = numpy.ndarray(3)    #: Cell basis vector 1
      self.cB1.fill(0)
      self.cB2 = numpy.ndarray(3)    #: Cell basis vector 2
      self.cB2.fill(0)
      self.cB3 = numpy.ndarray(3)    #: Cell basis vector 3
      self.cB3.fill(0)
      self.cO = numpy.ndarray(3)     #: Cell origin
      self.cO.fill(0)
      ############################################################################################

      ############################################################################################
      self.remang = 0                #: Remove angular momentum 
      self.remcom = 0                #: Remove Center of Mass motion?
      self.time = 0                  #: Current time in seconds
      self.temperature = 300         #: Kelvin temperature
      self.masses = numpy.ndarray(0)      #: Diagonal mass matrix (M)
      self.invmasses = numpy.ndarray(0)   #: Diagonal inverse mass matrix (M^-1)
      self.masssum = 0                    #: Sum over all atomic masses
      ###########################################################################################
      
      ###########################################################################################
      # accept and reject are used for Monte Carlo propagators
      # These are used to sample some metastable structures of the molecule but are non-physical
      # Given a current structure of the molecule, these pick another random structure and accept
      # or reject it based on a probability, usually e^-deltaE/kT where deltaE is the change in
      # energy between the old and new structure (the lower the energy change, the more likely 
      # the new structure is to be used otherwise we keep the old).
      # Both are false by default.
      self.accept = 0
      self.reject = 0
      ###########################################################################################
      
      ProtoMolApp.ProtoMolApp.turnOffHints() # Turn off hints in the SWIG-wrapped ProtoMol code
      # NOTE: self.positions and self.velocities are also
      # available numpy arrays.
      #####################################################################

      ############################################
      # NOTE: USER SHOULD NOT TOUCH THESE!
      # THESE ARE NUMPY ARRAY WRAPPERS
      # USER SHOULD ACCESS THROUGH
      # self.positions and self.velocities ONLY.
      # These classes are from the ProtoMol back end
      self.__dict__['myTop'] = GenericTopology.T_Periodic() # 'mytop' is a GenericTopology.T_Periodic()
      self.__dict__['posvec'] = Vector3DBlock.Vector3DBlock() # 'posvec' is  a Vector3DBlock.Vector3DBlock()
      self.__dict__['velvec'] = Vector3DBlock.Vector3DBlock() # 'velvec' is a Vector3DBlock.Vector3DBlock()
      self.__dict__['myPAR'] = PARReader.PAR() # 'myPAR' is a PARReader.PAR().  Data from a CHARMM parameter file (force field)
      self.__dict__['myPSF'] = PSFReader.PSF() # 'myPSF' is a PSFReader.PSF().  Data from a CHARMM Protein Structure File (what is bonded to what?)
      self.__dict__['myPDB'] = PDBReader.PDB() # 'myPDB' is a PDBReader.PDB()   Data from a Protein Data Bank file (atomic positions)
      self.__dict__['myXYZ'] = XYZReader.XYZ() # 'myXYZ' is a XYZReader.XYZ()   Data from an XYZ file (simple three-column format for
                                               # positions and/or initial velocities if they are not random.
      self.__dict__['myEig'] = EigenvectorReader.EigenvectorInfo() # myEig' is a EigenvectorReader.EigenvectorInfo(), used for Normal Modes.
      ############################################

      self.positions = numpy.ndarray(0)   #: Atomic position vector
      self.velocities = numpy.ndarray(0)  #: Atomic velocity vector

      self.dirty = 1   #: Dirty bit, signifies that the Physical structure has changed since the last build. 
                       #  If the dirty bit is set, we rebuild the Physical structure (call build()) before calculating
                       # any new forces
   
   # Code to run when a Physical object loses scope.
   # We do not automatically delete all its data members (which could be referred to elsewhere)
   # Thus we pass and 
   def __del__(self):
      pass 

   ###########################################################################
   # If we print a Physical object, this gets displayed to the screen
   def __str__(self):
      return "Physical.Physical"
   def __repr__(self):
      return "Physical.Physical"
   ###########################################################################

   # We want to have a way of copying MDL Physical objects without
   # setting internal variables equal to each other and then having changes
   # affect one another.
   # So we can now do phys2 = phys.copy()
   # This returns a brand new Physical object and copies every individual
   # variable between them, element-by-element
   def copy(self, retval, forces="", dt=-1):
      """
      Perform a deep copy, avoid reference assignment
      """
      retval.seed = self.seed # set seed to retval.seed
      retval.exclude = self.exclude # set exclude to retval.exclude
      retval.cellsize = self.cellsize # set cellsize to retval.cellsize
      retval.bc = self.bc # set bc to retval.bc
      retval.defaultCBV = self.defaultCBV # set defaultCBV to retval.defaultCBV
      retval.remcom = -1 # set -1 to retval.remcom
      retval.remang = -1 # set -1 to retval.remang
      retval.defaultCBV = self.defaultCBV #set defaultCBV to retval.defaultCBV
      retval.myTop = self.myTop # set myTop to retval.myTop
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

      if (dt != -1):
         retval.app = ProtoMolApp.ProtoMolApp()
         f = Forces.Forces()
         print "Making app" 
         retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, f.energies, dt)
         apps.append(retval.app)
      #retval.build()
      #return retval
    
   # As an alternative to copy(), you can do a deepcopy() and pass a different set of forces
   # and timestep.  This way you can copy the same physical system, but use a different
   # set of forces and dt value when running the simulation.  Note that we have to reassemble
   # the ProtoMolApp (SWIG-wrapped structure from ProtoMol). 
   def deepcopy(self, forces, dt):
         retval = self.copy()
         retval.app = ProtoMolApp.ProtoMolApp()
         print "Making app 2"
         retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, forces.energies, dt)

      
   # SPECIAL ACCESSOR FOR self.positions or self.velocities
   # TO GET DATA FROM WRAPPERS
   # phys.__getattr__ is called by default when you are getting a member variable of the Physical class.
   # i.e. when you take phys.somevariable Python calls phys.__getattr__("somevariable")
   # Internally, Python maps variable names to values in the __dict__ array, i.e.
   # __dict__["somevariable"] = somevalue.
   # But we have a special case with positions and velocities.
   # phys.positions has to return the internal array of positions, but as an MDVec.
   # Similar with phys.velocities.
   # getC() returns a C array of floating point values that is used in the ProtoMol code as well, but
   # we convert it to an MDVec for these purposes.
   def __getattr__(self, name):
      if (name == 'positions'): # if the name = 'positions'
        if (self.gpu == True):
          return self.gpuPositions
        else:
          return (self.__dict__['posvec'].getC()).view(MDVec) # return the position vector
      elif (name == 'velocities'):# if the name = 'velocities'
        if (self.gpu == True):
          return self.gpuVelocities
        else:
          return (self.__dict__['velvec'].getC()).view(MDVec)  # return the velocity vector
      elif (name == 'invmasses'):
        if (self.gpu == True):
          return self.gpuInvmasses
        else:
          return self.__dict__['self.invmasses']
      elif (name == 'time'): # if the name = 'time'
         return self.myTop.time # return myTop time
      else: # if name is not set to 'positions', 'velocities'  or 'time' then 
         print name # print name
         return self.__dict__[name] 

   # SPECIAL ASSIGNMENT FOR self.positions or self.velocities
   # TO SET DATA IN WRAPPERS
   # __setattr__ is called uatomatically when you set an attribute of something of type Physical
   # i.e. phys.seed = 1234 would actually do phys.__setattr__('seed',1234)
   # Note we have to know if this is the first time we have set an attribute of the class, if
   # that's the case we have to rebuild the Physical object.
   def __setattr__(self, name, val):
      ###################################################################################
      # Is this the first time setting an attribute?  If so we must build (see below)
      firsttime = False # firsttime is set to False
      if (not self.__dict__.has_key(name)): # if self.__dict__ does not have key(name) 
         firsttime = True # firsttime is set to True
      ###################################################################################
      # If we are setting phys.positions, we must do an element-by-element copy of the array
      # that it is being set to (stored in val)
      if (name == 'positions'): # if name = 'positions'
        if (self.gpu == True):
#          if (type(val)
          if (type(val) != numpy.ndarray):
#            print "Positions:\n",val.matrix
#            x=raw_input(' ')
#            print "VERY BAD THINGS - POSITIONS: ", type(val)
#            try:
#              print 1
#              print val, val.matrix
#            self.gpuPositions.set_matrix(val)
#            except AttributeError:
#              print 2
#              print val
#              self.gpuPositions.set_matrix(numpy.asarray(val, config.floatX))
#          else:
          #if (val != self):
            self.gpuPositions.set_matrix(val)
#          print "gpuPos: ", self.gpuPositions.matrix[0]
#          print "PosVec: ", self.__dict__['posvec'].getC()[0]
#          print "Type of gpuPositions.matrix: ", type(self.gpuPositions.matrix)
#          print "gpuPositions.matrix:\n", self.gpuPositions.matrix
#          print "Type of gpuPositions contents: ", type(self.gpuPositions.matrix[0])
#	  print "Calling setC for positions"
#          self.__dict__['posvec'].setC(self.gpuPositions.matrix, True)
	  
#          print (self.__dict__['posvec'].getC()[0] == self.gpuPositions.matrix[0])
#          print self.__dict__['posvec'].getC()[0], " : ", self.gpuPositions.matrix[0] 
#          print self.__dict__['posvec'].getC()[1], " : ", self.gpuPositions.matrix[1]
#          print self.__dict__['posvec'].getC()[2], " : ", self.gpuPositions.matrix[2]
#          print "Positions Matrix set!\n", self.__dict__['posvec'].getC()
#          raw_input()
        else:
          if type(val) != 'numpy.ndarray':
            self.__dict__['self.positions'] = val
          else:
            for i in range (0, len(self.positions)): # loop from 0 to the length of positions
              self.positions[i] = val[i] # set val[i] to positions[i]
            self.__dict__['posvec'].setC(val, False)
   
      # Similar case for the velocities.
      elif (name == 'velocities'): # if name = 'velocities'
        if (self.gpu == True):
          if (type(val) != numpy.ndarray):
#            print "Velocities:\n", val.matrix
#            raw_input()
#            try:
            #self.gpuVelocities.set_matrix(val)
#            except AttributeError:
#              self.gpuVelocities.set_matrix(val)
#          else:
#          if (type(val) != numpy.ndarray):
          #if (val != self):
            self.gpuVelocities.set_matrix(val)
#          print "gpuVel: ", self.gpuVelocities[0]
#          print "VelVec: ", self.__dict__['velvec'].getC()[0]
#	  print "Calling setC for velocities"
#          self.__dict__['velvec'].setC(self.gpuVelocities.matrix, True)
        else:
          if type(val) != 'numpy.ndarray':
            self.__dict__['self.velocities'] = val
          else:
            for i in range (0, len(self.__dict__['self.velocities'])): # loop from 0 to length of velocities
              self.velocities[i] = val[i] # set val[i] to velocities[i]
            self.__dict__['velvec'].setC(val, False)
      elif (name == 'gpu'):
        if (val == True):
          import pycuda.driver as cuda
#          from theano import shared
          from GPUMatrixPyC import GpuMatFunctions, GpuMatrix
          blanklist = numpy.zeros(self.numAtoms() * 3, numpy.float32)
#          g_a = cuda.mem_alloc(blanklist.nbytes)
#          g_b = cuda.mem_alloc(blanklist.nbytes)
          functions = GpuMatFunctions(blanklist, blanklist)
          p = numpy.asarray(self.__dict__['posvec'].getC(), numpy.float32)
          g_pos = cuda.mem_alloc(p.nbytes)
          cuda.memcpy_htod(g_pos, p)
          self.gpuPositions = GpuMatrix(g_pos, cuda.mem_alloc(p.nbytes), blanklist.copy(), functions, self.numAtoms() * 3)
          v = numpy.asarray(self.__dict__['velvec'].getC(), numpy.float32)
          g_vel = cuda.mem_alloc(v.nbytes)
          cuda.memcpy_htod(g_vel, v)
          self.gpuVelocities = GpuMatrix(g_vel, cuda.mem_alloc(v.nbytes), blanklist.copy(), functions, self.numAtoms() * 3)
          i = numpy.asarray(self.__dict__['self.invmasses'], numpy.float32)
          g_inv = cuda.mem_alloc(i.nbytes)
          cuda.memcpy_htod(g_inv, i)
          self.gpuInvmasses = GpuMatrix(g_inv, cuda.mem_alloc(i.nbytes), blanklist.copy(), functions, self.numAtoms() * 3)
#          self.gpuPositions = GpuMatrix(shared(self.__dict__['posvec'].getC()), functions)
#          self.gpuVelocities = GpuMatrix(shared(self.__dict__['velvec'].getC()), functions)
#          self.gpuInvmasses = GpuMatrix(shared(numpy.asarray(self.__dict__['self.invmasses'], numpy.float32)), functions)
        self.__dict__['gpu'] = val
      #similar case for inverse mass
      elif (name == 'invmasses'):
        if (self.gpu == True):
          if (type(val) != numpy.ndarray):
#            self.gpuInvmasses = val
#          else:
            self.gpuInvmasses.set_matrix(val)
        else:
          if type(val) != 'numpy.ndarray':
            self.__dict__['self.invmasses'] = val
          else:
            for i in range (0, len(self.invmasses)): # loop from 0 to length of velocities
              self.invmasses[i] = val[i] # set val[i] to invmass[i]

  #    elif (name == 'masses'):
#	  if type(val) != 'numpy.ndarray':
 #                self.__dict__['self.masses'] = val
  #        else:
   #              for i in range (0, len(self.masses)): # loop from 0 to length of velocities
    #                    self.masses[i] = val[i] # set val[i] to invmass[i]

	  
      # Here we are setting the boundary conditions
      # They will either be Vacuum or Periodic
      # Vacuum assumes no external forces, Periodic assumes wraparound and is generally
      # used for solvated systems.
      # The SWIG-wrapped structure for the topology of the system is T_Vacuum() for Vacuum
      # and T_Periodic() for periodic - so we have to set myTop to one of those
      elif (name == 'bc'): # if the name = 'bc'
         self.__dict__['bc'] = val         # then set val to self.__dict__['bc']
         if (val == "Vacuum"): # is val = "Vacuum" 
           self.myTop = GenericTopology.T_Vacuum() # then set GenericTopology.T_Vacuum() to myTop
         else: # if val is not a "Vacuum"
           self.myTop = GenericTopology.T_Periodic() # set GenericTopology.T_Periodic() to myTop
         if (not firsttime): # is this is not firsttime
            self.build() # build()
      # Real biologiyal simulation time
      elif (name == 'time'): # if name = 'time'
         val /= Constants.invTimeFactor() # val = val/~.0205
         self.myTop.time = val # set val to myTop time
      # These are all simple, we set their entry in the dictionary to the passed value and rebuild
      # if it's not the first time
      elif (name == 'seed' or  # is name is 'seed','exclude','cellsize','remcom', or 'remang'
            name == 'exclude' or
            name == 'cellsize' or
            name == 'remcom' or
            name == 'remang'):
         self.__dict__[name] = val # set val equal to the name
         if (not firsttime): # if not firsttime
            self.build() # build
      else: # any other name - set it's entry in the dictionary but do not rebuild - since the
	    # attribute is not a known member of the Physical class.
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
   # This is in units kcal/mol * A^-3.
   # We compute this using the following equation: 
   # P = (Exx + Eyy + Ezz) / 3V  *  (unit conversion)
   # Exx, Eyy and Ezz are based on the Virial Theorem:
   # http://en.wikipedia.org/wiki/Virial_theorem
   # This relates the average total kinetic energy with the
   # average total potential energy.
   # If we add this Virial term to the energy calculation (line 1), we
   # can then apply the ideal gas law:
   # PV = NRT
   # But since Average Kinetic Energy <KE> is equal to T/3R, then:
   # PV = N<KE>/3, and 
   # P = N<KE>/3V.  Then we compute <KE> as the sum of the three projections
   # of the virial term.
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
   # Computes this as a function of positions
   # You can do this as a function of atomic positions
   # Compute the bounding box (xmin, ymin, zmin) and (xmax, ymax, zmax) as the
   # min and max x,y,z over all atoms.  
   # Then V = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
   def volume(self):
      """
      Volume of the system.

      @rtype: float
      @return: System volume
      """
      return self.myTop.getVolume(self.posvec) # return system volume


   # SIZE OF THE SYSTEM
   # Simply calls numAtomms() function
   def N(self):
      """
      Number of atoms.

      @rtype: int
      @return: Number of atoms.
      """
      return self.numAtoms() # return Number of atoms
   
   ##############################################################################
   # The next set of functions returns data obtained from a CHARMM Protein
   # Structure File (PSF), which describes the topology of the molecule:
   # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node21.html
   # 
   # Note that this data never changes throughout a simulation.
   # This includes:
   # - Number of atoms (Function: numAtoms)
   # - Number of bonds (Function: numBonds)
   # - Number of angles (Function: numAngles)
   # - Number of dihedrals (Function: numDihedrals)
   # - Number of impropers (Function: numImpropers)
   # - Number of H+ donors (Function: numDonors)
   # - Number of H+ acceptors (Function: numAcceptors)
   ##############################################################################
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


   ##############################################################################
   # The next set of functions again returns data obtained from a CHARMM Protein
   # Structure File (PSF), which describes the topology of the molecule:
   # http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node21.html
   # 
   # Note that this data never changes throughout a simulation.
   # This set of functions, however - returns a specific element of the molecule.
   # This uses the classes at the top of the file, which describe these elements
   # and their attributes.
   # The classes are: Atom, Bond, Angles, Dihedral, Improper, Donor and Acceptor
   # Note that the value returned is a call to the appropriate constructor of the
   # respective class, and has the correct number of attributes.
   #
   # Each of these functions accepts one parameter for the index (starting from 1)
   # If I wanted information on atom 1, I would thus do phys.atom(1) - note that
   # we must convert that to an index starting from zero to work with the above classes
   # Thus we always use index-1 in the constructor.
   #
   # Functions:
   # - atom(index): get atom #index
   # - bond(index): get bond #index
   # - angle(index): get angle #index
   # - dihedral(index): get dihedral #index
   # - improper(index): get improper #index
   # - donor(index): get donor #index
   # - acceptor(index): get acceptor #index
   #
   # Also two others, for useful simulation data:
   # - mass(index): atomic mass for atom #index
   # - charge(index): atomic charge for atom #index
   ##############################################################################

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
   def bond(self, index):
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
   def angle(self, index):
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
   def dihedral(self, index):
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
   def improper(self, index):
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
   def donor(self, index):
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
   def acceptor(self, index):
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
      return self.myPSF.getAttributeReal("atom", atom-1, "charge") * Constants.sqrtCoulombConstant() # return the Atomic charge [e]
   #####################################################################################################
   
  
   # SYSTEM TEMPERATURE.
   # We compute thisusing the equation:
   # PV = NRT
   # Above, we computed P = N<KE>/3V
   # Thus: N<KE>/3 = NRT
   #        <KE> = 3RT
   #        T = <KE>/3R
   # To get the average kinetic energy, sum 0.5*m*v^2 for all atoms and divide by N (number of atoms)
   # That is why we need velocities to compute temperature.
   def temperature(self):
      """
      System temperature (K)

      @rtype: float
      @return: Kelvin temperature
      """
      return TopologyUtilities.temperature(self.myTop, self.velvec) # return the Kelvin temperature

   # Dihedral #index in radians
   # We need this for some computations, for instantce the Finite Temperature String Method
   # that recognizes a conformation of a molecule by its backbone dihedral angle values.
   def dihedralRad(self, index):
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

   # Initialize velocities to random values.
   # note this will depend on temperature, the higher the tmpereature the higher the possible range
   # of initial velocity values (since heat generally results in an increase in motion).
   def randomVelocity(self, T):
      """
      Assign random velocities.

      @type T: float
      @param T: Kelvin temperature.
      """
      TopologyUtilities.randomVelocity(T, self.myTop, self.velvec, self.seed) # get the randomVelocity of TopologyUtilities

   # This function updates the center of mass and anglular momentum of the molecule
   # Generally, this should be done after every update to positions and velocities, since the positions affect
   # the center of mass and the velocities affect angular momentum
   # Since this is an inconvenience for the user to have to worry about, we generally call this behind the
   # scenes whenever either vector is changed.
   def updateCOM_Momenta(self):
      """
      Update center of mass and angular momentum
      """
      TopologyUtilities.buildMolecularCenterOfMass(self.posvec,self.myTop) # update center of mass
      TopologyUtilities.buildMolecularMomentum(self.velvec,self.myTop) #update angular momentum

   # build() function for the Physical class
   # Its job is to update all "Behind the scenes" data that the user will not touch
   # Much of this data is SWIG-wrapped information from the ProtoMol code
   # Attribute changes in the Physical class will most often result in a call to this function
   # to initialize everything.
   def build(self):
      """
      Build the physical data.
      """
      ###################################################################################
      # TIME
      tm = -1 # tm is set to -1
      if (hasattr(self, "myTop")): # if self has the attribute "myTop"
          tm = self.myTop.time # tm = myTop.time
      ###################################################################################
      # BOUNDARY CONDITIONS
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
             i = 0
             while (i < numpy.size(self.positions)):
                if (self.positions[i] < v1[0]):
                   v1[0] = float(self.positions[i])
                if (self.positions[i] > v2[0]):
                   v2[0] = float(self.positions[i])
                if (self.positions[i+1] < v1[1]):
                   v1[1] = float(self.positions[i+1])
                if (self.positions[i+1] > v2[1]):
                   v2[1] = float(self.positions[i+1])
                if (self.positions[i+2] < v1[2]):
                   v1[2] = float(self.positions[i+2])
                if (self.positions[i+2] > v2[2]):
                   v2[2] = float(self.positions[i+2])
                i += 3
             v4.fill(Constants.periodicBoundaryTolerance()/2.0)
             v1 = v1 - v4
             v2 = v2 + v4
             v3 = v2 - v1
             self.cB1[0] = v3[0]
             self.cB2[1] = v3[1]
             self.cB3[2] = v3[2]
             self.cO = v1 + v3 * 0.5
             self.myTop.setBC(self.cB1[0],self.cB1[1],self.cB1[2],self.cB2[0],self.cB2[1],self.cB2[2],self.cB3[0],self.cB3[1],self.cB3[2],self.cO[0],self.cO[1],self.cO[2])
      ###################################################################################
      # STRUCTURE OF THE SYSTEM (PSF)
      self.myTop.setExclusion(self.exclude)
      #print self.myPSF.numAtoms()
      #print hasattr(self.myPAR, 'readFlag')
      if (self.myPSF.numAtoms() > 0 and hasattr(self.myPAR, 'readFlag')):
	 print "BUILDING TOPOLOGY..."
         GenericTopology.buildTopology(self.myTop, self.myPSF, self.myPAR, 0, self.myTop.makeSCPISM())
	 print "DONE BUILDING TOPOLOGY."
	 if (numpy.size(self.velocities) == 0):
            # This function actually returns a value, since it is
            # unused the runtime environment detects a memory leak.
            # disown() removes this concern.
            MathUtilities.randomNumber(self.seed) 
            # NOTE TO SELF: THIS WAS COMMENTED OUT FOR PM 3
            # IT MAY EFFECT RANDOM NUMBER CONSISTENCY
            #aaa = MathUtilities.randomNumberFirst(self.seed, 1)
	    print "USING RANDOM VELOCITIES..."
            TopologyUtilities.randomVelocity(self.temperature, self.myTop, self.velvec, self.seed)
	 if (self.remcom >= 0):
	   TopologyUtilities.removeLinearMomentum(self.velvec, self.myTop).disown()
         if (self.remang >= 0):
	   TopologyUtilities.removeAngularMomentum(self.posvec, self.velvec, self.myTop).disown()
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
      ###################################################################################

      ###################################################################################
      # SET SELFICAL TIME
      if (tm != -1): # if tm does not equal -1
          self.myTop.time = tm # set tm to myTop time
      ###################################################################################

      self.dirty = 0  # clean now!


