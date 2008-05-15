import sys
import BondForce
import AngleForce
import DihedralForce
#import DihedralIIIForce
import HarmDihedralForce
import ImproperForce
#import WrapperMetaForce
import SimpleFullForce
#import FullForce
import CutoffForce
#import FullEwaldForce
#import PMEwaldForce
#import MultiGridForce

    
class ForceFactory:
  def __init__(self):
    """
    Initializes mappings from boundary conditions (for bonded forces)
    and boundary conditions, algorithms and switching functions (for nonbonded
    forces) to SWIG-wrapped force object constructors (not instances, saving
    memory).
    """
    self.bondForces = {'Vacuum':BondForce.BSF_Vacuum,
                       'Periodic':BondForce.BSF_Periodic} #: Mapping from boundary conditions to two-atom bond force object constructors

    self.angleForces = {'Vacuum':AngleForce.ASF_VBC,
                        'Periodic':AngleForce.ASF_PBC} #: Mapping from boundary conditions to three-atom angle force object constructors

    self.dihedralForces = {'Vacuum':DihedralForce.DSF_Vacuum,
                           'Periodic':DihedralForce.DSF_Periodic} #: Mapping from boundary conditions to four-atom dihedral force object constructors

    self.improperForces = {'Vacuum':ImproperForce.ISF_Vacuum,
                           'Periodic':ImproperForce.ISF_Periodic} #: Mapping from boundary conditions to four-atom improper force object constructors

    self.harmDihedralForces = {'Vacuum':HarmDihedralForce.HDSF_Vacuum,
                               'Periodic':HarmDihedralForce.HDSF_Periodic} #: Mapping from boundary conditions to harmonic dihedral force object constructors

    self.hd = 0  #: Number of harmonic dihedral forces
    
    ######################################################################################################
    # THE FOLLOWING NONBONDED FORCE OBJECTS ARE PREDEFINED AND WRAPPED AS SHARED OBJECTS FOR EFFICIENCY  #
    # OBJECTS ARE REGISTERED IN THE FACTORY AS MULTI-LEVEL MAPPINGS                                      #
    # INDEX ONE: BOUNDARY CONDITIONS                                                                     #
    # INDEX TWO: ALGORITHM                                                                               #
    # INDEX THREE: SWITCHING FUNCTION(S)                                                                 #                                  
    ######################################################################################################
    self.ljForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_L},
                               'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_LJF,
                                         'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_LJF,
                                         'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_LJF,
                                         'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_LJF,
                                         'CmpCnCn':CutoffForce.NCSF_CCM_OAPVBC_CCNCNSF_LJF
                                         }
                               },
                     'Periodic':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_P_U_L},
                                 'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPPBC_C1SF_LJF,
                                           'C2':CutoffForce.NCSF_CCM_OAPPBC_C2SF_LJF,
                                           'Cutoff':CutoffForce.NCSF_CCM_OAPPBC_CSF_LJF,
                                           'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_LJF,
                                           'CmpCnCn':CutoffForce.NCSF_CCM_OAPPBC_CCNCNSF_LJF
                                           }
                                 }
                     }  #: Maps boundary conditions, algorithm and switching function to van der Waals force object constructor            

               
    self.cdeForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_CDE},
                                'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CFDE,
                                          'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CFDE,
                                          'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CFDE,
                                          'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CFDE,
                                          'ComplementCn':CutoffForce.NCSF_CCM_OAPVBC_CMPCNNSF_CFDE}},
                      'Periodic':{'Cutoff':{'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_CFDE}}
                      } #: Maps boundary conditions, algorithm and switching function to coulomb dielectric force object constructor - used for implicit solvation
                                               

    self.coulombForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_C},
                                    'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CF,
                                              'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CF,
                                              'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CF,
                                              'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CF,
                                              'ComplementCn':CutoffForce.NCSF_CCM_OAPVBC_CCNCNSF_CF},
                                    'SCPISM':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CSCPF,
                                              'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CSCPF,
                                              'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CSCPF,
                                              'CmpCnCn':CutoffForce.NCSF_CCM_OAPVBC_CCNSF_CSCPF,
                                              'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CSCPF}},
                          
                          'Periodic':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_P_U_C},
                                      'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPPBC_C1SF_CF,
                                                'C2':CutoffForce.NCSF_CCM_OAPPBC_C2SF_CF,
                                                'Cutoff':CutoffForce.NCSF_CCM_OAPPBC_CSF_CF,
                                                'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_CF,
                                                'ComplementCn':CutoffForce.NCSF_CCM_OAPPBC_CCNCNSF_CF}}
                                      

                          } #: Maps boundary conditions, algorithm and switching function to electrostatic force object constructor.  For fast electrostatics, an additional mapping is performed for the terms calculated (i.e. for Ewald, real reciprocal and correction).

    self.bornForces = {'Vacuum':{'Cutoff':{'C1':CutoffForce.NCBF_CCM_OAPVBC_C1SF_CBF,
                                           'C2':CutoffForce.NCBF_CCM_OAPVBC_C2SF_CBF,
                                           'Cn':CutoffForce.NCBF_CCM_OAPVBC_CNSF_CBF,
                                           'CmpCnCn':CutoffForce.NCBF_CCM_OAPVBC_CCNSF_CBF,
                                           'Cutoff':CutoffForce.NCBF_CCM_OAPVBC_CSF_CBF}}}

               
    self.ljCoulombForces = {'Vacuum':{'SimpleFull':{('Universal', 'Universal'):SimpleFullForce.NSFSF_V_U_L_U_C},
                                      'Cutoff':{('C2', 'C1'):CutoffForce.NCSF_CCM_OAPTVBC_C2SF_LJF_C1SF_CF,
                                                ('C2', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_C2SF_LJF_CSF_CF,
                                                ('Cn', 'Cn'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_CN_CF,
                                                ('Cn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_C_CF,
                                                ('CmpCnCn', 'CmpCnCn'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_CN_CF,
                                                ('CmpCnCn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_C_CF}

                                      },
                            'Periodic':{'SimpleFull':{('Universal', 'Universal'):SimpleFullForce.NSFSF_P_U_L_U_C},
                                        'Cutoff':{('C2', 'C1'):CutoffForce.NCSF_CCM_OAPTPBC_C2SF_LJF_C1SF_CF,
                                                  ('C2', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_C2SF_LJF_CSF_CF,
                                                  ('Cn', 'Cn'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_CN_CF,
                                                  ('Cn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_C_CF,
                                                  ('CmpCnCn', 'CmpCnCn'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_CN_CF,
                                                  ('CmpCnCn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_C_CF}

                                        }
                            } #: Maps boundary conditions, algorithm and switching function pair to a unifed van der Waals and electrostatic force object constructor.  This saves performance by determining atom pairs just once for both types of pairwise forces.


    #########################################
    # WRAPPERS
    #########################################
    #self.mollify = WrapperMetaForce.WrapperMetaForce #: Wrapper function for MOLLY forces.  Wraps either a bond or angle force object.

   
  def createBondForce(self,bc):
    """
    Return a bond force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped bond force object.
    """
    return self.bondForces[bc]()
      
  def createAngleForce(self, bc):
    """
    Return an angle force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped angle force object.
    """
    return self.angleForces[bc]()
   
  def createDihedralForce(self, bc):
    """
    Return a dihedral force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped dihedral force object.
    """
    return self.dihedralForces[bc]()

  def createImproperForce(self, bc):
    """
    Return an improper force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped improper force object.
    """
    return self.improperForces[bc]()

  def createHarmDihedralForce(self, bc, params):
    """
    Return a harmonic dihedral force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @type params: dict
    @param params: Mapping from parameter names to corresponding values.

    @rtype: Force
    @return: SWIG-wrapped harmonic dihedral force object.
    """
    if (str(type(params['kbias']))[7:11] == 'list'):
      self.hd += 1  # INCREASE HARMONIC DIHEDRAL COUNTER
      return self.harmDihedralForces[bc](params['kbias'][self.hd-1], params['dihedralnum'][self.hd-1], params['angle'][self.hd-1], 0)
    else:
      return self.harmDihedralForces[bc](params['kbias'], params['dihedralnum'], params['angle'], 0)      

  def createMollyBondForce(self, bc):
    """
    Return a mollified bond force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped mollified bond force object.
    """
    return self.mollify("MollyBond",1,self.bondForces[bc](),"MollyBond")

  def createMollyAngleForce(self, bc):
    """
    Return a mollified angle force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped mollified angle force object.
    """
    return self.mollify("MollyAngle",1,self.angleForces[bc](),"MollyAngle")

  def createLennardJonesForce(self, bc, params):
      """
      Return a van der Waals force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped van der Waals force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal")        
      newforce = self.ljForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params)

  def createCoulombDiElecForce(self, bc, params):
      """
      Return an coulomb dielectric force object (for implicit solvation).
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped coulomb dielectric force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal") 
      newforce = self.cdeForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params)
       	
  def createCoulombForce(self, bc, params, fastelectro):
      """
      Return an electrostatic force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped electrostatic force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal")
         
      if (alg == "Ewald"):
         real = self.getParameter(fastelectro['Ewald'], 'real', True)
         reciprocal = self.getParameter(fastelectro['Ewald'], 'reciprocal', True)
         correction = self.getParameter(fastelectro['Ewald'], 'correction', True)
         terms = str(real)[0].lower()+str(reciprocal)[0].lower()+str(correction)[0].lower()
         newforce = self.coulombForces[bc][alg][switch][terms]()
      elif (alg == "PME"):
         real = self.getParameter(fastelectro['PME'], 'real', True)
         reciprocal = self.getParameter(fastelectro['PME'], 'reciprocal', True)
         correction = self.getParameter(fastelectro['PME'], 'correction', True)
         interp = self.getParameter(fastelectro['PME'], 'interpolation', "BSpline")
         terms = str(real)[0].lower()+str(reciprocal)[0].lower()+str(correction)[0].lower()         
         newforce = self.coulombForces[bc][alg][switch][interp][terms]()
      elif (alg == "MultiGrid"):
         direct = self.getParameter(fastelectro['MultiGrid'], 'direct', True)
         correction = self.getParameter(fastelectro['MultiGrid'], 'correction', True)
         smooth = self.getParameter(fastelectro['MultiGrid'], 'smooth', True)
         interp = self.getParameter(fastelectro['MultiGrid'], 'interpolation', "BSpline")
         terms = str(direct)[0].lower()+str(correction)[0].lower()+str(smooth)[0].lower()
         newforce = self.coulombForces[bc][alg][switch][interp][terms]()
      else:
        newforce = self.coulombForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params, fastelectro)

  def createBornForce(self, bc, params):
      alg = "Cutoff"
      switch = self.getParameter(params, 'switching', "Universal")
      newforce = self.bornForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params)
  #def createLennardJonesTableForce(self, bc, params):
  #    alg = "Cutoff"  # Must be the case
  #    switch = self.switchingFunctions(params, 2)
  #    newforce = ljTableForces[bc][alg][switch]
  #    if (params.has_key('cutoff') and params['cutoff'] != -1):
  #             newforce.setCutoff(params['cutoff'])

  def createLennardJonesCoulombForce(self, bc, params):
      """
      Return a coupled van der Waals - electrostatic force evaluation object.
      Saves performance by calculating atom pairs only once for both
      type of forces.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped coupled van der Waals and electrostatic force object.
      """
      alg = self.getParameter(params, 'algorithm', "Cutoff")  # Different default
      switch = self.switchingFunctions(params, 2)
      newforce = self.ljCoulombForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params)
      
  #def createLennardJonesCoulombERForce(self, bc, params):
  #    alg = "Cutoff"  # Must be the case
  #    switch = self.switchingFunctions(params, 2)        
  #    newforce = self.ljCoulombERForces[bc][alg][switch]()      
  #    return self.applyParameters(newforce, bc, alg, switch, params)

      
  #def createLennardJonesCoulombERTForce(self, bc, params):
  #    alg = "Cutoff"  # Must be the case
  #    switch = self.switchingFunctions(params, 2)   
  #    newforce = self.ljCoulombERTForces[bc][alg][switch]()
  #    return self.applyParameters(newforce, bc, alg, switch, params)
  
  #def createLennardJonesCoulombMGDForce(self, bc, params):
  #    alg = "Cutoff"  # Different default here
  #    switch = self.switchingFunctions(params, 3) 
  #    newforce = self.ljCoulombMGDForces[bc][alg][switch]()
  #    return self.applyParameters(newforce, bc, alg, switch, params)
       
  #def createLennardJonesCoulombMGDTForce(self, bc, params):
  #    alg = "Cutoff"  # Must be the case
  #    switch = self.switchingFunctions(params, 3) 
  #    newforce = self.ljCoulombMGDTForces[bc][alg][switch]()
  #    return self.applyParameters(newforce, bc, alg, switch, params)

  def createMagneticDipoleForce(self, bc):
      """
      Return a magnetic dipole force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)

      @rtype: Force
      @return: SWIG-wrapped magnetic dipole force object.
      """
      newforce = self.magneticDipoleForces[bc]()
      return newforce


  def switchingFunctions(self, params, number):
      """
      Little helper function, to take a set of parameters
      and obtain a specific number of switching functions,
      setting defaults as necessary (i.e. if no switching functions
      are provided, Universal is assumed)

      @type params: dict
      @param params: Mapping from parameter names to values

      @type number: int
      @param number: Number of switching functions to obtain
      """
      fxns = ()
      if (not params.has_key('switching')):
        for i in range(0, number):
           fxns += ("Universal",)
      else:
        if (number == 1):
           fxns += (params['switching'],)
        else:
           if (str(type(params['switching']))[7:11] == 'list'):
             for j in range(0, len(params['switching'])):
                fxns += (params['switching'][j],)
             for k in range(len(params['switching']), number):
                fxns += (fxns[k-1],)
           else:
             for j in range(0, number):
               fxns += (params['switching'],)  # Just add the first
      return fxns

  # Return the parameter if it exists, otherwise return the provided default value
  def getParameter(self, params, name, defaultval=0):
    """
    Return a parameter value if it exists, otherwise return the passed default value.

    @type params: dict
    @param params: Mapping from parameter names to values

    @type name: string
    @param name: Name of the parameter

    @type defaultval: (any type)
    @param defaultval: Default value for the parameter if it is not found
    """
    if (params.has_key(name)):
      return params[name]
    else:
      return defaultval
    
  # Add application of parameters here
  # This is its own function, because a lot of algorithm-switching function combinations
  # will have the same extra parameters, and this avoids code repetition
  def applyParameters(self, newforce, bc, alg, switch, params, fastelectro=None):
    """
    Depending on the algorithm used, set parameter values for the passed
    force object.

    @type newforce: Force
    @param newforce: Pairwise force object

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @type alg: string
    @param alg: Algorithm for pairwise evaluation

    @type switch: string
    @param switch: Switching function

    @type params: dict
    @param params: Mapping from parameter names to values

    @type fastelectro: dict
    @param fastelectro: Mapping from fast electrostatic parameter names to values

    @rtype: Force
    @return: Pairwise force object with instantiated parameters.
    """

    if (str(newforce).find('CoulombForceDiElec') != -1): # Dielectric forces have extra params
      eps = self.getParameter(params, 'epsilon', 1)
      d = self.getParameter(params, 'd', 78)
      s = self.getParameter(params, 's', 0.3)
      if (alg == "SimpleFull"):
        newforce = newforce.makeNewDiElec(self.getParameter(params, 'blocksize', 32), eps, d, s)
      elif (alg == "Cutoff"):
        if (switch == "C1"):
          newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'), eps, d, s)
        elif (switch == "C2"):
          newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'),
                                            eps, d, s,
                                            self.getParameter(params, 'switchon'))
        else:
          newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'),
                                            eps, d, s,
                                            self.getParameter(params, 'switchon'),
                                            self.getParameter(params, 'switchoff'),
                                            self.getParameter(params, 'order', 2))      
    else:
      if (alg == "SimpleFull"):
        newforce = newforce.makeNew(self.getParameter(params, 'blocksize', 32))
      elif (alg == "Cutoff"):
        if (switch == "C1" or switch == "Cutoff"):
          newforce = newforce.makeNew(self.getParameter(params, 'cutoff'))
        elif (switch == "C2"):
          newforce = newforce.makeNew(self.getParameter(params, 'cutoff'),
                                      self.getParameter(params, 'switchon'))
        else:
          newforce = newforce.makeNew(self.getParameter(params, 'cutoff'),
                                      self.getParameter(params, 'switchon'),
                                      self.getParameter(params, 'switchoff'),
                                      self.getParameter(params, 'order', 2))
    return newforce
