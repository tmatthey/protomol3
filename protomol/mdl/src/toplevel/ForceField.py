from ForceGroup import *
import sys

class ForceField(ForceGroup):
   # DEFAULT DESTRUCTOR WILL DELETE FORCE ARRAY
   # THIS DESTRUCTION IS HANDLED IN THE BACK END
   # SO WE DON'T DO IT AGAIN
   def __del__(self):
      """
      Destructor.
      """
      return

   def setDefaults(self):
      """
      Set default values for all parameters, for all force evaluation algorithms.  This includes bonded and nonbonded, pairwise and fast electrostatic.
      """
      self.params = {

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

       'MagneticDipole':{'algorithm':'Cutoff',
                         'switching':'C1',
                         'cutoff':0},

       # Harmonic Dihedral forces are an exception, in that
       # there can be more than one per force evaluation.
       # i.e., the user may wish to constrain different dihedrals.
       # Thus each parameter is a list of values.
       'HarmonicDihedral':{'kbias':[],          # scaling factor
                           'dihedralnum':[],    # dihedral index
                           'angle':[]}          # constraint angle in radians
       } #: Mapping from parameter names to default values for bonded and nonbonded forces, except fast electrostatics



      # Only used if the user chooses a fast electrostatic algorithm
      # for Coulombic forces.
      self.fastelectro = {'Ewald':{'real':True,                     # real term
                                   'reciprocal':True,               # reciprocal term
                                   'correction':True,               # correction term
                                   'alpha':'optimal',               # splitting parameter
                                   'accuracy':0.00001,              # accuracy
                                   'expansion':3.0},                # expansion factor

                          'PME':{'real':True,                       # real term
                                 'reciprocal':True,                 # reciprocal term
                                 'correction':True,                 # correction term
                                 'alpha':'optimal',                 # splitting parameter
                                 'accuracy':0.00001,                # accuracy
                                 'cutoff':6.0,
                                 'expansion':3.0,                   # expansion factor
                                 'interpolation':'BSpline',         # BSpline or Hermite
                                 'gridsize':5},                     # grid dimensions
                          
                          'MultiGrid':{'direct':True,               # direct term
                                       'correction':True,           # correction term
                                       'smooth':True,               # smoothing term
                                       'interpolation':'Hermite',   # BSpline or Hermite
                                       'levels':2,                  # inter/anterpolation levels
                                       's':10,                      # softening distance
                                       'order':6,                   # interpolation order
                                       'ratio':2,                   # fine-to-coarse ratio
                                       'coarsegridsize':8,          # size of toplevel grid
                                       'finegridsize':3,            # size of finest grid (VBC)
                                       'finegridorigin':[0,0,0]}    # origin of finest grid (VBC)
                          } #: Mapping from parameter names to default values for fast electrostatics


      # A dirty bit.  If this is set, on a call to propagate()
      # the forces will be rebuilt with the above parameters.
      # Upon each rebuild this bit is cleared, upon each changing
      # of any of the above parameters the bit is set.
      self.dirty = 1  #: Dirty bit, allows for lazy building of force objects upon system propagation

         
      ##################################################
      # USER-ACCESSIBLE STRUCTURES

      ###################################################################
      # ARRAY OF FORCE TYPES.
      # CAN CONTAIN THE FOLLOWING VALUES:
      # b - BOND, a - ANGLE, d - DIHEDRAL, i - IMPROPER
      # l - LENNARDJONES, c - COULOMB
      # lc - LENNARDJONES/COULOMB TOGETHER PAIRWISE EVALUATION (Default)
      # m - MAGNETICDIPOLE
      # mb - MOLLY BOND, ma - MOLLY ANGLE
      self.forcetypes = [] #: List of force types, can contain 'b' (bond), 'a' (angle), 'd' (dihedral), 'i' (improper), 'h' (harmonic dihedral), 'l' (van der Waals), 'c' (electrostatic), 'e' (implicit solvation dielectric scaling), 'lc' (coupled vdW and electrostatic), 'm' (magnetic dipole), 'mb' (mollified bond) and 'ma' (mollified angle)
      self.pythonforces = []
      ################################################################### 
      
      ###################################################################
      # DEFAULT WILL CREATE AN EMPTY FORCE FIELD, BUT IF THE CHARMM
      # FLAG IS SET, WE MUST ADD THE SIX CORE FORCES
      if (self.charmm):
         self.bondedForces("badi")
         self.nonbondedForces("lc")
      ###################################################################   

   def __setattr__(self, s, t):
      if (s == 'params' or s == 'fastelectro'):   # IF WE CHANGE params or fastelectro, it's dirty
         self.dirty = 1
      self.__dict__[s] = t
      
   # REMOVE ALL BONDED FORCES FROM THE EVALUATION
   def removeBondedForces(self):
      """
      Remove all bonded (bond, angle, dihedral, improper, harmonic dihedral)
      forces from the force field.
      """
      pos = self.findForce('b')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('a')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('d')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('i')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = 0
      while (pos != -1):
         pos = self.findForce('h')
         if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)

   # REMOVE ALL NONBONDED FORCES FROM THE EVALUATION
   def removeNonbondedForces(self):
      """
      Remove all nonbonded forces (van der Waals, electrostatic, magnetic dipole) from the force field.
      """
      pos = self.findForce('l')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('c')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('lc')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('m')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('e')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)


   # REMOVE ALL MOLLY FORCES FROM THE EVALUATION
   def removeMollyForces(self):
      """
      Remove all mollified forces from the force field.
      """
      pos = self.findForce('mb')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)
      pos = self.findForce('ma')
      if (pos != -1):
          self.forcetypes.__delslice__(pos, pos+1)

      
   # ADD BONDED FORCES ACCORDING TO THE PASSED STRING
   # b=BOND, a=ANGLE, d=DIHEDRAL, i=IMPROPER
   # EXAMPLE INVOCATION: bondedForces("ba")
   def bondedForces(self, inputstring):
      """
      Add bonded forces contained in the input string ('b', 'a', 'd', 'i', or 'h')

      @type inputstring: string
      @param inputstring: Contains the characters representing the bonded forces to instantiate
      """
      self.removeBondedForces()
      if (inputstring.find('b') != -1):
          self.forcetypes.append('b')
      if (inputstring.find('a') != -1):
          self.forcetypes.append('a')
      if (inputstring.find('d') != -1):
          self.forcetypes.append('d')
      if (inputstring.find('i') != -1):
          self.forcetypes.append('i')
      if (inputstring.count('h') != 0):
          for i in range(0, inputstring.count('h')):
            self.forcetypes.append('h')
      #SystemFactory.buildForces(self.theforces)


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
          inputstring.find('c') != -1):
          self.forcetypes.append('lc')
      else:
          if (inputstring.find('l') != -1):
              self.forcetypes.append('l')
          if (inputstring.find('c') != -1):
              self.forcetypes.append('c')
      if (inputstring.find('m') != -1):
          self.forcetypes.append('m')
      if (inputstring.find('e') != -1):
          self.forcetypes.append('e')
      #SystemFactory.buildForces(self.theforces)


   # ADD MOLLY FORCES ACCORDING TO THE PASSED STRING
   # b=BOND, a=ANGLE, d=DIHEDRAL, i=IMPROPER
   # EXAMPLE INVOCATION: mollyForces("ba")  
   def mollyForces(self, inputstring):
      """
      Add mollified forces contained in the input string, ('b' or 'a') are available.

      @type inputstring: string
      @param inputstring: Contains the characters representing the mollified forces to instantiate
      """
      self.removeMollyForces()
      if (inputstring.find('b') != -1):
          self.forcetypes.append('mb')
      if (inputstring.find('a') != -1):
          self.forcetypes.append('ma')


   # FIND A FORCE AND RETURN ITS INDEX
   # IF NOT FOUND, RETURN -1.
   def findForce(self, t):
      """
      Search for a force in the array of force types.

      @type t: char
      @param t: Type of force ('b', 'a', etc. see above).  Force can be bonded, nonbonded, or molly.

      @rtype: int
      @return: Index of this type of force in the types array; -1 if not found.
      """
      ii = 0
      while (ii < self.forcetypes.__len__()):
          if (self.forcetypes[ii] == t):
              return ii
          ii += 1
      return -1


   def addPythonForce(self, pyforce):
      self.pythonforces.append(pyforce)
      #self.addSystemForce(PySystemForce(pyforce))

   # BREAK THE LENNARDJONES/COULOMB EVALUATION OF THE PAIRS
   # TO SEPARATE COMPUTATIONS.  IF THEY ARE ALREADY SEPARATE
   # DO NOTHING
   def breakLennardJonesCoulombForce(self):
      """
      Breaks a coupled van der Waals - electrostatic force into individual calculations.  This is invoked if the user specifies a different algorithm for van der Waals and electrostatic - for example if one is direct and the other uses a cutoff; they will not use the same set of atom pairs.
      """
      pos = self.findForce('lc')
      if (pos != -1):
         self.forcetypes.remove(self.forcetypes[pos])
         self.forcetypes.append('l')
         self.forcetypes.append('c')

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
      if (forces.forcevec.size() != phys.numAtoms()):
         forces.forcevec.resize(phys.numAtoms())
      forces.forcevec.zero()
      forces.energies.clear()
      self.evaluateSystemForces(phys.app, forces.forcevec)
      self.evaluateExtendedForces(phys.app, forces.forcevec)
   def calculateMollyForces(self, phys, forces):
      """
      Similar to calculateForces, but for mollified forces.

      @type phys: Physical
      @param phys: The Physical system.

      @type forces: Forces
      @param forces: MDL Forces object.
      """
      self.evaluateMollyForces(phys.myTop, phys.positions, forces.force, forces.energies)

