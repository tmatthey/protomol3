import MathUtilities
from PropagatorFactory import *
import time
import Constants
import TopologyUtilities
import numpy
import ProtoMolApp
import copy
import ModifierShake

class Propagator:
   """
   Provides functionality to propagate a system with time.
   """
   def __init__(self, phys, forces, io):
      #####################################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.myPropagator = 0                  #: Propagator object
      self.myStep = 0                        #: Current simulation step
      self.myTimestep = 0                    #: Propagator timestep
      self.myLevel = 0                       #: Number of levels
      #####################################################################################

      phys.build() #Build the physical data

      self.phys = phys #: Physical object 
      self.forces = forces #: Forces object
      self.io = io  #: IO object
      io.phys = phys

      
   def reset(self):
      """
      Reset the state of the Propagator object.
      """
      self.myPropagator = 0 # Propagator object set to zero
      self.myStep = 0  # current simulation step to zero
      self.myTimestep = 0 # Propagator timestep to zero
      self.myLevel = 0 # Number of levels to zero

   def isMDL(self, integ):
      """
      Determine whether or not a propagator has been coded in MDL.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @rtype: boolean
      @return: True if the passed propagator is prototyped in pure Python.
      """
      return (hasattr(integ, "prerunmodifiers"))# check for the existence of integ propagator having prerunmodifiers


   def addPreInitModifier(self, integ, modifier):
      """
      Add a modifier to execute before propagator initialization.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.preinitmodifiers.append(modifier)# append modifier functions to preinitmodifiers (before propagator initialization)


   def addPostInitModifier(self, integ, modifier):
      """
      Add a modifier to execute after propagator initialization.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postinitmodifiers.append(modifier) # append modifier functions to postinitmodifiers (after propagator initialization)


   def addPreRunModifier(self, integ, modifier):
      """
      Add a modifier to execute before propagator execution.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.prerunmodifiers.append(modifier) # append modifier functions to prerunmodifiers (before propagator execution)


   def addPostRunModifier(self, integ, modifier):
      """
      Add a modifier to execute after propagator execution.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postrunmodifiers.append(modifier) # append modifier functions to postrunmodifiers (after propagator execution)



   def addPreForceModifier(self, integ, modifier):
      """
      Add a modifier to execute before force calculation.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.preforcemodifiers.append(modifier) # append modifier to preforcemodifiers (before force calculation)


   def addPostForceModifier(self, integ, modifier):
      """
      Add a modifier to execute after force calculation.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postforcemodifiers.append(modifier) # append modifier to postforcemodifiers (after force calculation)


   # RUN A LIST OF PASSED MODIFIERS ON THE PASSED PROPAGATOR
   def runModifiers(self, modifiers, phys, forces, prop, integ):
      """
      Run modifiers of a propagator

      @type modifiers: list of functions
      @param modifiers: A set of routines which alternates propagator behavior

      @type phys: Physical
      @param phys: MDL Physical object

      @type forces: Forces
      @param forces: MDL Forces object

      @type prop: Propagator
      @param prop: MDL Propagator object
      
      @type integ: object
      @param integ: MDL propagator object (STS or MTS)
      """
      #integ.postf
      for ii in range(0, modifiers.__len__()): # for-loop(0,length of modifiers)
         modifiers[ii](phys, forces, prop, integ) #loop through modifiers array to run each individual modifier of propagation

   def timestep(self, integ):
      """
      Return the timestep of a propagator, scaled accordingly

      @type integ: object
      @param integ: MDL propagator object (STS or MTS)

      @rtype: float
      @return: The timestep (dt) of a propagator
      """
      return integ.getTimestep() * Constants.invTimeFactor() # multiply timestep by ~.0205

   def calculateForces(self, forces):
        """
        Calculate forces and update the atomic force vector.
      
        @type forces: Forces
        @param forces: MDL Forces object
        """
        for ii in range(0, self.myPropagator.preforcemodifiers.__len__()): # for-loop(0, length of preforcemodifiers)
           self.myPropagator.preforcemodifiers[ii](self.phys, forces, self, self.myPropagator) # run each preforcemodifier
        self.myPropagator.calculateForces() # calculate forces of myPropagator
        forces.forcevec = self.myPropagator.getForces()# set myPropagator equivalent to forevec
        for ii in range(0, self.myPropagator.postforcemodifiers.__len__()):# for-loop(0, length of postforcemodifiers)
           self.myPropagator.postforcemodifiers[ii](self.phys, forces, self, self.myPropagator) # run each postforcemodifier 

   def initNext(self, phys, forces):
      """
      For multiple timestepping, initialize the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      """
      tempI = self.myPropagator # set myPropagator to temp variable
      self.myLevel += 1 # increase myLevel by +=1
      setPropagator(self, phys, forces, self.myPropagator.next, levelswitch=True) # run setPropagator
      self.myPropagator = tempI # set temp variable to myPropagator
      self.myLevel -= 1 # decrease myLevel by +=1

   def runNext(self, phys, forces, cL):
      """
      For multiple timestepping, execute the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      @type cL: int
      @param cL: Cycle length (number of times to execute the
                 inner propagator)
      """
      tempI = self.myPropagator # set myPropagator to temp variable
      self.myPropagator = self.myPropagator.next # set myPropagator to my Propagator.next
      self.myLevel = self.myLevel + 1 # set myLevel to myLevel + 1
      executePropagator(self, phys, forces, self.io, cL) # run executePropagator
      self.myLevel = self.myLevel - 1 # decrease myLevel by 1
      self.myPropagator = tempI # set temp variable to myPropagator


   def finishNext(self, phys, forces, prop):
      """
      For multiple timestepping, finish the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      @type prop: Propagator
      @param prop: MDL Propagator object

      """
      tempI = self.myPropagator # set myPropagator to temp variable
      self.myPropagator = self.myPropagator.next # set myPropagator to myPropagator.next
      if (self.isMDL(self.myPropagator)): # if myPropagator isMDL
         self.myPropagator.finish(phys, forces, prop) #finish myPropagator
      self.myPropagator = tempI # set temp variable to myPropagator

   def deepCopyForce(self, ff): # ff is forcefield
      p = ff.params # set forcefield parameters to p
      if (ff.charmm): # if forcefield.charmm
         ff2 = self.forces.makeForceField(self.phys, "charmm") # set makeForceField(charmm) to ff2
      else: # else if not charmm
         ff2 = self.forces.makeForceField(self.phys) # set makeForceField(physical) to ff2
      for ii in ff.forcetypes: # loop through each element in the forcetype array 
         ff2.forcetypes.append(ii) # append ii to forcetypes
      ff2.params = p # set p to ff2.params
      ff2.build() # build ff2
      #self.forces.removeForceField(ff)
      import ForceGroup
      ForceGroup._swig_setattr_nondynamic(ff2, ForceGroup.ForceGroup, "thisown", 0) # this function keeps ff2 and ff from 
                                                                                    # resetting back to original defaults 
      return ff2 # return ff2

   # PROPAGATE THE SYSTEM
   # USE METHOD "name"
   # arg1 = NUMBER OF STEPS
   # arg2 = TIMESTEP
   # arg3 = ForceField STRUCTURE
   # args = OPTIONAL EXTRA ARGUMENTS AS TUPLES
   def propagate(self, scheme="Leapfrog", steps=0, cyclelength=-1, dt=0.1, forcefield=[], params={}):
       """
       Propagate the system.

       @type scheme: string
       @param scheme: Name of the propagator to use.
       
       @type steps: int
       @param steps: Number of steps for execution.

       @type cyclelength: int
       @param cyclelength: Cyclelength for MTS propagation (-1 is STS).
       
       @type dt: float
       @param dt: Timestep.

       @type forcefield: ForceField
       @param forcefield: MDL ForceField object.

       @type params: dictionary
       @param params: Extra parameters unique for this propagation scheme.
                     (This could be empty).

       """
       self.myTimestep = dt # set myTimestep to dt
       chain = ()
       if (cyclelength != -1):  # MTS
          if (str(type(cyclelength))[7:11] == 'list'): # LIST, MANY LEVELS
             levels = len(cyclelength) + 1 # levels are equal to the length of the cyclelength +1
             outertime = cyclelength[0] # outertime is equal to the first element in cyclelength
          else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
             levels = 2 # levels set to 2
             outertime = cyclelength # outertime set to cyclelength

          if (str(type(scheme))[7:11] == 'list'): # LIST, MANY LEVELS
             outerscheme = scheme[0] # outerscheme equal to the first element in scheme
          else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
             outerscheme = scheme # outerscheme equal to scheme


          # THE NUMBER OF FORCEFIELDS PROVIDED MUST EQUAL THE NUMBER
          # OF PROPAGATION LEVELS
          if (len(forcefield) != levels): # if the length of forcefield is not equivalent to levels
             print "[MDL] Error in propagate(): ", levels, " levels of propagation with ", len(forcefield), " force fields."
          outerforcefield = self.deepCopyForce(forcefield[0]) # print this out

          if (str(type(scheme))[7:11] != 'list'): # if the type of scheme is not a list MTS
             chain += (params,) # chain = chain +(params,) 
          else: # if a list
             if (params.has_key(outerscheme)): # True if params has a key outerscheme
                 chain += (params[outerscheme],)# chain = chain + params[outerscheme]
             else: # if params does not have a key outerscheme
                 chain += ({},) # chain = chain +({},)
	  for i in range(1, levels): # for i in range from 1, levels
             if (str(type(scheme))[7:11] == 'list' and i < len(scheme)): # if the type of scheme is a list and i is less than the
									 # length of scheme	
                chain += (scheme[i],) # chain = chain + (scheme[i],)
             if (str(type(cyclelength))[7:11] == 'list' and i < len(cyclelength)): # if the type of cyclelength is a list and
								            # i is less than the length of cyclelength
                chain += (cyclelength[i],) # chain = chain + (cyclelength[i],)
             else: # if not a list
                chain += (dt,) # chain = chain + (dt,)
	     chain += (self.deepCopyForce(forcefield[i]),) #chain = chain +(deepCopyForce(forcefield[i]),)
             if params.has_key(scheme[i]): # if params has key scheme[i] 
                chain += (params[scheme[i]],) # add chain to (params[scheme[i]],)
             else: # if params has no key scheme[i]
                chain += ({},) # chain = chain +({},)
       else: #STS
          outertime = dt # set the outertime to dt
          outerscheme = scheme # set outerscheme to scheme
          outerforcefield = self.deepCopyForce(forcefield) # set outerforcefield to deepCopyForce(forcefield)
          chain += (params,) # chain = chain + (params,)
       # Build force fields.
       # Tricky, because we could be dealing with
       # a single object or list.
       if ((str(type(forcefield)))[7:len(str(type(forcefield)))-2] == 'list'): # if the type of forcefield is a list (MTS)
          for ff in forcefield: # step through each element in forcefield
             if (ff.dirty): ff.build() # and if ff has a dirty bit then ff.build()
	     if (ff.gbsa): # if ff.gbsa
	        self.phys.myTop.implicitSolvent = 2 # set implicitSolvent to 2
		self.phys.myTop.doGBSAOpenMM = 1 # set doGBSAOpenMM to 1
		self.phys.build() # build
       else: # if the type of forcefield is not a list 
          if (forcefield.dirty): # if forcefield has a dirty bit
              forcefield.build() # then build forcefield
	  if (forcefield.gbsa): # if forcefield.gbsa
	      self.phys.myTop.implicitSolvent = 2 # set implicitSolvent to 2
              self.phys.myTop.doGBSAOpenMM = 1 # set doGBSAOpenMM to 1
	      self.phys.build() # physical build()
       if (self.io.dirty): # if self.io is using the dirty bit
          self.io.build() # io.build()
        

       if (propFactory.getType(outerscheme) == "method"): # if the type fo outerscheme is equal to a "method"
          # Calculate the forces, store them in force.
          if (not hasattr(self.phys, 'app')): # if self.phys does not have the attribute 'app'
             self.phys.app = ProtoMolApp.ProtoMolApp() # set phys.app to the ProtomolApp()
             self.phys.app.makeApp(self.phys.myTop, self.phys.posvec, self.phys.velvec, self.forces.energies, dt) # makeApp
          self.phys.updateCOM_Momenta() # update the center of mass and angular momentum
          outerforcefield.calculateForces(self.phys, self.forces) # calculate the forces in outerforcefield
          self.io.run(self.phys, self.forces, 0, outertime) # run 
          self.io.myProp = self # set self to myProp
          for ii in range(1, steps+1): # for-loop ii, (1,steps+1)
             self.phys.app.energies.clear() # clear the energies in phys.app
             self.forces.energies.initialize(self.phys) # initialize the forces with self.phys energies
             propFactory.create(1, outerscheme, self.phys, self.forces, self.io, 1, outertime*Constants.invTimeFactor(), outerforcefield, *chain) # create
             self.phys.time = ii*outertime*Constants.invTimeFactor() # set the physical time to loop variable* outertime* ~.0205 
             self.io.run(self.phys, self.forces, ii, outertime) # run
             self.phys.updateCOM_Momenta() # update Center of mass and angular momentum
             #self.phys.app.energies.clear()
             #self.forces.forcevec.zero()
       else: # Object
	  setPropagator(self, self.phys, self.forces, propFactory.applyModifiers(propFactory.create(1, outerscheme, outertime, outerforcefield, *chain), outerscheme)) # set propagator
          shake = False # shake is false
          if (params.has_key('shake') and params['shake'] == 'on'): # if params has the key 'shake' and params['shake'] is turned on
              shake = True # shake is set to true
              shakeMod = ModifierShake.ModifierShake(0.000001, 30) # .000001 is set to epsilon and max iteration is 30
              #shakeMod = self.myPropagator.createShakeModifier(0.000001, 30)
              self.myPropagator.adoptPostDriftOrNextModifier(shakeMod) #set shakeMod to adoptPostDriftOrNextModifier
          rattle = False # rattle is false
          if (params.has_key('rattle') and params['rattle'] == 'on'): # if params has key 'rattle' and params['rattle'] is turned on
              rattle = True # set rattle to true
              rattleMod = self.myPropagator.createRattleModifier(0.02, 30) # epsilon is .02 and max iteration is set to 30
              self.myPropagator.adoptPostDriftOrNextModifier(rattleMod) #set rattleMod to adoptPostDriftOrNextModifier
          executePropagator(self, self.phys, self.forces, self.io, steps)  # Runs the propagator for a number of steps
          if (shake): # if shake
             self.myPropagator.removeModifier(shakeMod) # remove shakemod in myPropagator
          if (rattle): # if rattle
             self.myPropagator.removeModifier(rattleMod)# remove rattlemod in myPropagator
       self.phys.updateCOM_Momenta() # update Center of Mass and angular momentum
           
