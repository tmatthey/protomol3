from Energies import *
import Vector3DBlock
from ForceField import *
from ForceFactory import *
import MathUtilities
import numpy
import Forces
from PySystemForce import *

class Forces:
   """
   Contains the atomic force vector, and a structure
   to hold system energies.
   """
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
      self.__dict__['forcevec'] = Vector3DBlock.Vector3DBlock() # set forcevec as a Vector3DBlock.Vector3DBlock()
      #self.force = 0  #: Atomic force vector.
      #########################################################
 
      self.accept = 0 # set accept and reject to 0
      self.reject = 0
         
   # SPECIAL ACCESSOR FOR self.forces
   # TO GET DATA FROM WRAPPERS   
   def __getattr__(self, name): # get attribute
      if (name == 'force'): # if the name is 'force'
         return self.__dict__['forcevec'].getC() # return the data of forcevec
      else: # or if not 'forces'
         return self.__dict__[name] # return other data

   # SPECIAL ASSIGNMENT FOR self.forces
   # TO SET DATA IN WRAPPERS 
   def __setattr__(self, name, val):# set attribute
      if (name == 'force' and self.__dict__.has_key('force')): # if the name is 'force' and the dict has 'force'
         self.__dict__['forcevec'].setC(val) # then set forcevec the specified value
      if type(val) != 'numpy.ndarray' and name == 'force':
	 self.__dict__['force'] = val	 
      else:
         self.__dict__[name] = val # set other data to another name

   def __str__(self):
      return "Forces.Forces"
   
   def __repr__(self):
      return "Forces.Forces"


   def initializeEnergies(self, app):
      """
      Initialize the energies structure.
      """
      self.energies = app.energies 
      
   # RESET SYSTEM STATE TO DEFAULTS
   def reset(self):
      """
      Reset data members to default values{myForcefields, energies and forces}.
      """
      self.myForceFields = list() # reset myForcefields by setting it to an empty list()
      self.energies = Energies() # set energies to Energies() which will reset energies
      self.forces.fill(0) # fill forces with 0 , resetting forces

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
      # this function simply sets up a new force field and then returns the forcefield initialized with 
      # forces and energies
      ff = ForceField() # set up ff as an empty ForceField(). 
      ff.charmm = 0 # set the charmm forces to 0
      ff.forceFactory = ForceFactory() # set up ff.forceFactory as a ForceFactory()
      if (args.__len__() > 0): # if the length of the arguments are greater than 0
         ff.charmm = 1 # set the charmm forces to 1
      ff.theforces = self # set self equivalent to "theforces"
      self.phys = ff.phys = phys # set phys
      ff.bc = phys.bc # set bc
      ff.setDefaults()# set the defaults for ff
      self.energies.initialize(phys)# now initialize energies with phys
      ff.thisown = 0 
      #self.myForceFields.append(ff)
      return ff # return the force field
      #return self.myForceFields[self.myForceFields.__len__()-1]


   # NOTE: THIS WILL PERMANENTLY DELETE THE PASSED FORCE FIELD OBJECT
   # FROM THE ARRAY!  USE WITH CAUTION.
   def removeForceField(self, ff):
      """
      Remove the passed force field from the member list.

      @type ff: ForceField
      @param ff: MDL force field.      
      """
      for i in range(0, self.myForceFields.__len__()): # goes from 0 to the length of myForceFields
         if (self.myForceFields[i] == ff):             # setting the forcefields at index i to ff
	    self.myForceFields.remove(ff)              # then deleting each individual forcefield at the specified index
            break                                      

   def randomForce(self, phys, seed):
    # we use random force for Langevin functions which use an N->V->T system. (N(number of Atoms, V(Volume) and T(temperature))
    # randomForce is then used in heat baths. (boiling water or heated solution) Langevin functions are a type of second order 
    # differential equation that uses a damping force. 
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

        ff = 0 # starts the index at 0        
        while (ff < phys.numAtoms()*3): # to make a random forcefield filled with random numbers
           retval[ff] = MathUtilities.randomGaussianNumber(seed)
           retval[ff+1] = MathUtilities.randomGaussianNumber(seed)
           retval[ff+2] = MathUtilities.randomGaussianNumber(seed)
           ff += 3
        return retval 

