from ForceGroup import *
import sys
import PySystemForce
class ForceField(ForceGroup):
   """
   A holder for the following data:
     1. Set of forces to evaluate.
     2. Algorithms to use for evaluation.
     3. Switching function(s) to apply, if any.
     4. Extra parameters, depending on the particular forces.
   """
   # DEFAULT DESTRUCTOR WILL DELETE FORCE ARRAY
   # THIS DESTRUCTION IS HANDLED IN THE BACK END
   # SO WE DON'T DO IT AGAIN
   #def __del__(self):
   #   """
   #   Destructor.
   #   """
   #   print "DEF DEL"
   #   pass

   def setDefaults(self):
      """
      Set default values for all parameters, for all force evaluation algorithms.  This includes bonded and nonbonded, pairwise and fast electrostatic.
      """
      self.params = { # set the parameters to include the LennardJones, Coulomb, LennardJonesCoulomb, CoulombDiElec and HarmonicDihedral forces

       'LennardJones':{'algorithm':'SimpleFull',
                       'switching':'Universal',
                       'cutoff':-1},

       'Coulomb':{'algorithm':'SimpleFull',
                  'switching':'Universal',
                  'cutoff':-1},

       # The 'type' parameter of 'LennardJonesCoulomb' can be set to:
       # 'original', 'EwaldReal', 'EwaldRealTable', 'MagneticDipole',
       # 'MagneticDipoleRealTable'
       'LennardJonesCoulomb':{'algorithm':'SimpleFull',
                              'switching':['Universal','Universal'],
                              'cutoff':-1,
                              'type':'original'},

       'CoulombDiElec':{'algorithm':'SimpleFull',
                        'switching':'Universal',
                        'cutoff':-1},

       # Harmonic Dihedral forces are an exception, in that
       # there can be more than one per force evaluation.
       # i.e., the user may wish to constrain different dihedrals.
       # Thus each parameter is a list of values.
       'HarmonicDihedral':{'kbias':[],          # scaling factor
                           'dihedralnum':[],    # dihedral index
                           'angle':[]}          # constraint angle in radians
       } #: Mapping from parameter names to default values for bonded and nonbonded forces, except fast electrostatics


      # A dirty bit.  If this is set, on a call to propagate()
      # the forces will be rebuilt with the above parameters.
      # Upon each rebuild this bit is cleared, upon each changing
      # of any of the above parameters the bit is set.
      self.dirty = 1  #: Dirty bit, allows for lazy building of force objects upon system propagation
      self.gbsa = False #: Flag telling whether or not to do generalized Born
         
      ##################################################
      # USER-ACCESSIBLE STRUCTURES

      ###################################################################
      # ARRAY OF FORCE TYPES.
      # CAN CONTAIN THE FOLLOWING VALUES:
      # b - BOND, a - ANGLE, d - DIHEDRAL, i - IMPROPER
      # l - LENNARDJONES, c - COULOMB
      # lc - LENNARDJONES/COULOMB TOGETHER PAIRWISE EVALUATION (Default)
      self.forcetypes = [] #: List of force types, can contain 'b' (bond), 'a' (angle), 'd' (dihedral), 'i' (improper), 'h' (harmonic dihedral), 'l' (van der Waals), 'c' (electrostatic), 'e' (implicit solvation dielectric scaling), 'lc' (coupled vdW and electrostatic)
      self.pythonforces = [] #: List of Python-prototyped orce objects
      ################################################################### 
      
      ###################################################################
      # DEFAULT WILL CREATE AN EMPTY FORCE FIELD, BUT IF THE CHARMM
      # FLAG IS SET, WE MUST ADD THE SIX CORE FORCES
      if (self.charmm):
         self.bondedForces("badi") # bonded forces: bond, angle, dihedral, and improper
         self.nonbondedForces("lc")# non-bonded forces: LennardJones and Coulomb 
      ###################################################################   

   def __setattr__(self, s, t):
      if (s == 'params'):   # IF WE CHANGE params or fastelectro, it's dirty
         self.dirty = 1 
      self.__dict__[s] = t 
      
   # REMOVE ALL BONDED FORCES FROM THE EVALUATION
   def removeBondedForces(self):
      """
      Remove all bonded (bond, angle, dihedral, improper, harmonic dihedral)
      forces from the force field.
      """
      # by setting each Bonded force to !=-1 it will reset to original defaults (in this case removing them)
      if (self.findForce('b') != -1): self.forcetypes.remove('b')
      if (self.findForce('a') != -1): self.forcetypes.remove('a')
      if (self.findForce('d') != -1): self.forcetypes.remove('d')
      if (self.findForce('i') != -1): self.forcetypes.remove('i')
      if (self.findForce('h') != -1):  # Can be multiple harmonic dihedral forces
         while (self.forcetypes.count('h') != 0): self.forcetypes.remove('h')


   # REMOVE ALL NONBONDED FORCES FROM THE EVALUATION
   def removeNonbondedForces(self):
      """
      Remove all nonbonded forces (van der Waals, electrostatic, magnetic dipole) from the force field.
      """
      # same scenario here by setting the Nonbonded forece !=-1 it will reset to original defaults
      if (self.findForce('l') != -1): self.forcetypes.remove('l')
      if (self.findForce('c') != -1): self.forcetypes.remove('c')
      if (self.findForce('lc') != -1): self.forcetypes.remove('lc')
      if (self.findForce('e') != -1): self.forcetypes.remove('e')

      
   # ADD BONDED FORCES ACCORDING TO THE PASSED STRING
   # b=BOND, a=ANGLE, d=DIHEDRAL, i=IMPROPER
   # EXAMPLE INVOCATION: bondedForces("ba")
   def bondedForces(self, inputstring):
      """
      Add bonded forces contained in the input string ('b', 'a', 'd', 'i', or 'h') 

      @type inputstring: string
      @param inputstring: Contains the characters representing the bonded forces to instantiate
      """
      self.removeBondedForces() # remove all the bonded forces 
      if (inputstring.find('b') != -1): self.forcetypes.append('b') # append bond force
      if (inputstring.find('a') != -1): self.forcetypes.append('a') # append angular force
      if (inputstring.find('d') != -1): self.forcetypes.append('d') # append dihedral forece
      if (inputstring.find('i') != -1): self.forcetypes.append('i') # append improper force
      if (inputstring.count('h') != 0): # there ared different harmonic bond forces 
          for i in range(0, inputstring.count('h')): self.forcetypes.append('h') # append all harmonic bond forces
      self.dirty = 1 # and set dirty to 1

   # ADD NONBONDED FORCES ACCORDING TO THE PASSED STRING
   # l=LENNARDJONES, c=COULOMB, m=MAGNETIC DIPOLE
   # EXAMPLE INVOCATION: nonbondedForces("lc")
   def nonbondedForces(self, inputstring):
      """
      Add nonbonded forces contained in the input string ('l', 'c', 'm', or 'e').  If 'l' and 'c' are included a coupled vdW-electrostatic is instantiated for maximum performance.

      @type inputstring: string
      @param inputstring: Contains the characters representing the nonbonded forces to instantiate
      """
      self.removeNonbondedForces()
      if (inputstring.find('l') != -1 and
          inputstring.find('c') != -1): self.forcetypes.append('lc') # appending the lennardcoulomb force
      else: # if the lennard and coulomb forces are not set both to -1
          if (inputstring.find('l') != -1): self.forcetypes.append('l') # append the lennard force 
          if (inputstring.find('c') != -1): self.forcetypes.append('c') # append the coulomb force
      if (inputstring.find('e') != -1): self.forcetypes.append('e') # append the electrostatic force
      self.dirty = 1 # and set dirty to 1

   # FIND A FORCE AND RETURN ITS INDEX
   # IF NOT FOUND, RETURN -1.
   def findForce(self, t):
      """
      Search for a force in the array of force types.

      @type t: char
      @param t: Type of force ('b', 'a', etc. see above).  Force can be bonded or nonbonded.

      @rtype: int
      @return: Index of this type of force in the types array; -1 if not found.
      """
      for ii in range(0, len(self.forcetypes)): # search from 0 to length of forcetypes for a specified force
          if (self.forcetypes[ii] == t): return ii # return the index of the Force
      return -1 # this means that the force was not found


   def addPythonForce(self, pyforce):
      """
      Add a Python-prototyped force object for evaluation.

      @type pyforce: PySystemForce
      @param pyforce: An instance of the Python-prototyped force.
      """
      self.pythonforces.append(pyforce) # append python force


   # BREAK THE LENNARDJONES/COULOMB EVALUATION OF THE PAIRS
   # TO SEPARATE COMPUTATIONS.  IF THEY ARE ALREADY SEPARATE
   # DO NOTHING
   def breakLennardJonesCoulombForce(self):
      """
      Breaks a coupled van der Waals - electrostatic force into individual calculations.  This is invoked if the user specifies a different algorithm for van der Waals and electrostatic - for example if one is direct and the other uses a cutoff; they will not use the same set of atom pairs.
      """
      pos = self.findForce('lc') # set pos to find the lennardcoulomb forece
      if (pos != -1): # if they are found
         self.forcetypes.remove('lc') # remove lennardcoulomb force
         self.forcetypes.append('l')  # and append the forces separatly lennard
         self.forcetypes.append('c')  # and coulomb force, this would break them apart.

   # THESE LAST TWO FUNCTIONS
   # CAN EVALUATE FORCES DIRECTLY THROUGH THE FORCEFIELD OBJECT
   # ONLY NEED TO USE THEM WHEN CONSTRUCTING A PROPAGATOR FUNCTION
   # (i.e. NOT A CLASS WHICH INHERITS FROM STS OR MTS).  OTHERWISE,
   # DO NOT EVEN USE THEM.
   def calculateForces(self, phys, forces):
      """
      Calculate all forces (except molly) in the force field, and populate the passed MDL forces object.  This should be called from a Python-prototyped propagation function.
      
      @type phys: Physical
      @param phys: The Physical system.

      @type forces: Forces
      @param forces: MDL Forces object.
      """
     # self.__dict__['forcevec'] = Vector3DBlock.Vector3DBlock()
      if (forces.forcevec.size() != phys.numAtoms()): # if the size of forcevec is not equivalent to the size of numAtoms
         forces.forcevec.resize(phys.numAtoms()) # resize 
      forces.forcevec.zero() # and set the forcevec of forces to zero
      #print len(phys.positions)
      #print phys.app.positions.size()
      phys.app.energies.clear() # clear energies
      self.phys = forces.phys = phys 
      #phys.posvec.setC(phys.positions)
      self.evaluateSystemForces(phys.app, forces.forcevec) # now evaluate the System forces: 
      #sys.exit(1)
      self.evaluateExtendedForces(phys.app, forces.forcevec)# and evaluate the Extended forces :

   def build(self):
    """
    Using an MDL force factory, instantiate all force objects stored in the forcetypes data member of each force field.  These will be SWIG-wrapped objects and are appended to the forcearray data member of each force field.
    """
    # TMC 1-13-08: I think this can be improved.  Future:
    # 1. Have a build() member of the ForceField class.
    # 2. Make the ForceFactory a singleton.
    # 3. Give the ForceFactory a method which maps force characters to creation functions.  In this way there are multiple mapping levels.
 
    self.forceFactory.hd = 0 # set the harmonic-dihedral force to 0

    if (self.params['LennardJones'] != # if the parameters of the LennardJones and Coulomb forces are not equivalent break them.
        self.params['Coulomb']):       # and set dirty to 0
            self.breakLennardJonesCoulombForce()
    self.dirty = 0

    self.forcearray = [] # reset forcearray

    bornflag=-1
    # the for-loop is used primarily to step through each element of forcetypes
    for forcetype in self.forcetypes:
          if (forcetype == 'b'): # if there exists a bond force in forcetype
              self.forcearray.append(self.forceFactory.createBondForce(self.bc)) # append and create a Bond force in forcearray
          elif (forcetype == 'a'): # if there exists an angular force in forcetype
              self.forcearray.append(self.forceFactory.createAngleForce(self.bc)) # append and create an Angle force in forcearray
          elif (forcetype == 'd'): # if there exists a dihedral force  in forcetype
              self.forcearray.append(self.forceFactory.createDihedralForce(self.bc)) # append and create a Dihedral Force in forcearray
          elif (forcetype == 'i'): # if there exists an improper force in forcetype
              self.forcearray.append(self.forceFactory.createImproperForce(self.bc)) # append and create an Improper Force in forcearray
          elif (forcetype == 'l'): # if there exists a lennard force in forcetype
              self.forcearray.append(self.forceFactory.createLennardJonesForce(self.bc, self.params['LennardJones'])) # append
                                                                                  # the lennard Jones force with the same parameters
          elif (forcetype == 'c'): # if there exists a coulomb force
              if (self.params['Coulomb']['algorithm'] == 'SCPISM'): # and if the parameters are as follows equivalent to SCPISM
                  # this means that it's in an "implicit Solvent"--> H20
                 self.phys.myTop.doSCPISM = 1 # set doSCPISM to 1 -> True, o ->False
                 self.phys.build() # build()
                 #if (not self.params['Coulomb'].has_key('NoBorn')):
                 #   print "CREATING BORN FORCE"
                 #   self.forcearray.append(self.forceFactory.createBornForce(self.bc, self.params['Coulomb']))
                 #   bornflag = len(self.forcearray)-1
                 #else:
                 #    print "NOT CREATING BORN FORCE"
              if (self.params['Coulomb']['algorithm'] == 'GB' or  # if the params  = Generalized Born or GBACE
	          self.params['Coulomb']['algorithm'] == 'GBACE'):
		  self.forcearray.append(self.forceFactory.createBornBurial()) # append and create BornBurial to forcearray
                  self.addForce(self.forcearray[len(self.forcearray)-1])# add force at the specified index of forcearray - 1
	          self.forcearray.append(self.forceFactory.createBornRadii()) # append and create BornRadii to forcearray 
                  self.addForce(self.forcearray[len(self.forcearray)-1])# add force at the specified index of forcearray - 1
              if (not self.params['Coulomb'].has_key('OnlyBorn')): # if the parameters of the Coulomb fore does not have 'OnlyBorn'      
                 self.forcearray.append(self.forceFactory.createCoulombForce(self.bc, self.params['Coulomb']))# append the coulomb force to forcearray 
          elif (forcetype == 'e'): # if the forcetype is electrostatic
              self.forcearray.append(self.forceFactory.createCoulombDiElecForce(self.bc, self.params['CoulombDiElec']))
              # append the DiElectric force to forcearray
          elif (forcetype == 'lc'): # if the forcetype is lennardJonesCoulomb
             self.forcearray.append(self.forceFactory.createLennardJonesCoulombForce(self.bc, self.params['LennardJonesCoulomb']))
              # append the LennardJonesCoulomb force to forcearray
          elif (forcetype == 'h'): # if the forcetype is harmonic 
              self.forcearray.append(self.forceFactory.createHarmDihedralForce(self.bc, self.params['HarmonicDihedral']))
              # append the Harmonic Dihedral force to forcearray
          self.addForce(self.forcearray[len(self.forcearray)-1])# add force at the specified index of forcearray - 1
	  import ForceGroup
	  ForceGroup._swig_setattr_nondynamic(self.forcearray[len(self.forcearray)-1], ForceGroup.Force, "thisown", 0) # this is 
                                                                                        # used to make sure the variables do not reset

    if (bornflag != -1): self.forcetypes.insert(bornflag, 'c') # if the bornflag is not at its default 
                                                               # forcetype[bornflag] = 'Coulomb'
    for pyforce in self.pythonforces: # iterate through pythonforces
         self.forcearray.append(PySystemForce.PySystemForce(pyforce)) #append pySystemForce(pyforce) to forcearray
	 PySystemForce._swig_setattr_nondynamic(self.forcearray[len(self.forcearray)-1], PySystemForce.PySystemForce, "thisown", 0)
         # used to make sure variable do not reset
         self.addSystemForce(self.forcearray[len(self.forcearray)-1]) # then finally addSystemForce at the end of forcearray -1
