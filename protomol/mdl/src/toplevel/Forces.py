from Energies import *
import Vector3DBlock
from ForceField import *
from ForceFactory import *
import MathUtilities
import numpy
from PySystemForce import *

class Forces:
   def __init__(self):
      #####################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.myForceFields = list()   #: Array of MDL force fields
      self.energies = Energies()    #: Holds system energies
      # NOTE: self.forces is a numpy array which holds system forces.
      #####################################################################

      #########################################################
      # NOTE: USER SHOULD NOT TOUCH THIS!
      # THIS IS A WRAPPER FOR THE FORCE ARRAY
      # WHICH IS ACCESSIBLE THROUGH self.force
      self.__dict__['forcevec'] = Vector3DBlock.Vector3DBlock()
      #self.force = 0  #: Atomic force vector.
      #########################################################

      
   def dirty(self):
      """
      Checks if any of the force fields created is dirty.

      @rtype: boolean
      @return: True if any of the force fields has been modified since the last propagation; false otherwise.
      """
      for ff in self.myForceFields:
         if (ff.dirty == 1):
            return True
      return False
         
   # SPECIAL ACCESSOR FOR self.forces
   # TO GET DATA FROM WRAPPERS   
   def __getattr__(self, name):
      if (name == 'force'):
         return self.__dict__['forcevec'].getC()
      else:
         return self.__dict__[name]

   # SPECIAL ASSIGNMENT FOR self.forces
   # TO SET DATA IN WRAPPERS 
   def __setattr__(self, name, val):
      if (name == 'force' and self.__dict__.has_key('force')):
         self.__dict__['forcevec'].setData(val)
      else:
         self.__dict__[name] = val

   # RESET SYSTEM STATE TO DEFAULTS
   def reset(self):
      """
      Reset data members to default values.
      """
      self.myForceFields = list()
      self.energies = Energies()
      self.forces.fill(0)

   # CREATE A FORCE FIELD AND APPEND IT TO THE END
   # OF THE ARRAY.
   # ALSO RETURN THE NEW FORCE FIELD POINTER.
   def makeForceField(self, phys, *args):
      """
      Create a new MDL force field.

      @type phys: Physical
      @param phys: The physical system.

      @type args: tuple
      @param args: Python tuple, may have no values or a string \"charmm\"

      @rtype: ForceField
      @return: Newly instantiated MDL force field.
      """
      ff = ForceField()
      ff.charmm = 0
      ff.forceFactory = ForceFactory()
      if (args.__len__() > 0):
         ff.charmm = 1
      ff.theforces = self
      self.phys = phys
      ff.bc = phys.bc
      ff.setDefaults()
      self.myForceFields.append(ff)
      return self.myForceFields[self.myForceFields.__len__()-1]


   # NOTE: THIS WILL PERMANENTLY DELETE THE PASSED FORCE FIELD OBJECT
   # FROM THE ARRAY!  USE WITH CAUTION.
   def removeForceField(self, ff):
      """
      Remove the passed force field from the member list.

      @type ff: ForceField
      @param ff: MDL force field.      
      """
      for i in range(0, self.myForceFields.__len__()):
         if (self.myForceFields[i] == ff):
            self.myForceFields.erase(ff)


   def randomForce(self, phys, seed):
        """
        Compute a random (x, y, z) force.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type seed: integer
        @param seed: Random number generator seed.

        @rtype: numpy.ndarray
        @return: The random force as a three-element array (x,y,z)
        """
        retval = numpy.ndarray(phys.numAtoms()*3)
        #tmp = numpy.ndarray(3)

        ff = 0        
        while (ff < phys.numAtoms()*3):
           retval[ff] = mathutilities.randomGaussianNumber(seed)
           retval[ff+1] = mathutilities.randomGaussianNumber(seed)
           retval[ff+2] = mathutilities.randomGaussianNumber(seed)
           #retval[ff:ff+3] = tmp*numpy.sqrt(phys.mass(ff/3+1))
           ff += 3
        return retval


   def build(self):
    """
    Using an MDL force factory, instantiate all force objects stored in the forcetypes data member of each force field.  These will be SWIG-wrapped objects and are appended to the forcearray data member of each force field.
    """
    # TMC 1-13-08: I think this can be improved.  Future:
    # 1. Have a build() member of the ForceField class.
    # 2. Make the ForceFactory a singleton.
    # 3. Give the ForceFactory a method which maps force characters to creation functions.  In this way there are multiple mapping levels.

 
    for forcefield in self.myForceFields:
      #forcefield.clear()
      forcefield.forceFactory.hd = 0

      if (forcefield.params['LennardJones'] !=
          forcefield.params['Coulomb']):
            forcefield.breakLennardJonesCoulombForce()
      forcefield.dirty = 0

      forcefield.forcearray = []

      for forcetype in forcefield.forcetypes:
          if (forcetype == 'b'):
              forcefield.forcearray.append(forcefield.forceFactory.createBondForce(forcefield.bc))
          elif (forcetype == 'a'):
              forcefield.forcearray.append(forcefield.forceFactory.createAngleForce(forcefield.bc))
          elif (forcetype == 'd'):
              forcefield.forcearray.append(forcefield.forceFactory.createDihedralForce(forcefield.bc))
          elif (forcetype == 'i'):
              forcefield.forcearray.append(forcefield.forceFactory.createImproperForce(forcefield.bc))
          elif (forcetype == 'l'):
              forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesForce(forcefield.bc, forcefield.params['LennardJones']))
          elif (forcetype == 'c'):
              if (forcefield.params['Coulomb']['algorithm'] == 'SCPISM'):
                 if (forcefield.params['Coulomb'].has_key('bornswitch')):
                    self.phys.myTop.doSCPISM = forcefield.params['Coulomb']['bornswitch']
                 else:
                    self.phys.myTop.doSCPISM = 1
                 self.phys.build()
                 if (not forcefield.params['Coulomb'].has_key('NoBorn')):
                    forcefield.forcearray.append(forcefield.forceFactory.createBornForce(forcefield.bc, forcefield.params['Coulomb']))
              if (not forcefield.params['Coulomb'].has_key('OnlyBorn')):      
                 forcefield.forcearray.append(forcefield.forceFactory.createCoulombForce(forcefield.bc, forcefield.params['Coulomb'], forcefield.fastelectro))
          elif (forcetype == 'e'):
              forcefield.forcearray.append(forcefield.forceFactory.createCoulombDiElecForce(forcefield.bc, forcefield.params['CoulombDiElec']))
          elif (forcetype == 'lc'):
              if (not forcefield.params['LennardJonesCoulomb'].has_key('type') or
                  forcefield.params['LennardJonesCoulomb']['type'] == 'original'):
                  forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesCoulombForce(forcefield.bc, forcefield.params['LennardJonesCoulomb']))
              elif (forcefield.params['LennardJonesCoulomb']['type'] == 'EwaldReal'):
                  forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesCoulombERForce(forcefield.bc, forcefield.params['LennardJonesCoulomb']))
              elif (forcefield.params['LennardJonesCoulomb']['type'] == 'EwaldRealTable'):
                  forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesCoulombERTForce(forcefield.bc, forcefield.params['LennardJonesCoulomb']))
              elif (forcefield.params['LennardJonesCoulomb']['type'] == 'MagneticDipole'):
                  forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesCoulombMGDForce(forcefield.bc, forcefield.params['LennardJonesCoulomb']))
              elif (forcefield.params['LennardJonesCoulomb']['type'] == 'MagneticDipoleTable'):
                  forcefield.forcearray.append(forcefield.forceFactory.createLennardJonesCoulombMGDTForce(forcefield.bc, forcefield.params['LennardJonesCoulomb']))
          elif (forcetype == 'mb'):
              forcefield.forcearray.append(forcefield.forceFactory.createMollyBondForce(forcefield.bc))
          elif (forcetype == 'ma'):
              forcefield.forcearray.append(forcefield.forceFactory.createMollyAngleForce(forcefield.bc))
          elif (forcetype == 'h'):
              forcefield.forcearray.append(forcefield.forceFactory.createHarmDihedralForce(forcefield.bc, forcefield.params['HarmonicDihedral']))
          elif (forcetype == 'm'):
              forcefield.forcearray.append(forcefield.forceFactory.createMagneticDipoleForce(forcefield.bc))

          # ADD TO THE BACK END
          if (forcetype == 'ma' or
              forcetype == 'mb'): # MOLLY FORCE
             forcefield.addMetaForce(forcefield.forcearray[forcefield.forcearray.__len__()-1])
          else: # SYSTEM FORCE
             forcefield.addForce(forcefield.forcearray[forcefield.forcearray.__len__()-1])
      for pyforce in forcefield.pythonforces:
         forcefield.forcearray.append(PySystemForce(pyforce))
         forcefield.addSystemForce(forcefield.forcearray[forcefield.forcearray.__len__()-1])
