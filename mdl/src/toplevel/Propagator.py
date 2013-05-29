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
       
       ##############################################################################################
       # This block of code assembles a propagator chain
       # In the case of multiple time-stepping, we could have an outer integrator with an arbitrary
       # set of arguments, and an inner integrator with a different set of arguments
       # The chain theoretically has no limit on how long it can be, although chains of length
       # larger than three are pretty rare
       chain = () # chain is initially empty
       if (cyclelength != -1):  # MTS if the cyclelength is NOT -1 (its default value)
	  # Now, it cyclelength is not -1, there are two possibilities:
	  # 1. cyclelength is a list of lengths (i.e., [3,2]) which means there is more than
	  #    one outer integrator
	  # 2. cyclelength is an integer (i.e., 3) in which case there is only one outer integrator
	  # Both cases have an inner integrator
          if (str(type(cyclelength))[7:11] == 'list'): # Case 1
             levels = len(cyclelength) + 1 # levels: The number of integrators we have: # of outers (length of list) plus the one inner integrator
             outertime = cyclelength[0] # outertime is equal to the cyclelength of the outermost integrator
          else: # Case 2, one outer integrator and one integrator
             levels = 2 # levels set to 2 (one outer, one inner)
             outertime = cyclelength # outertime set to cyclelength - there is only one outer integrator and that is the outermost

	  # Similar to the cyclelength case above, but this deals with schemes
	  # Again, this can be a list (i.e., ["Impulse, "Leapfrog"]) or just a string (i.e., "Leapfrog")
	  # We set outerscheme as the outermost scheme, which is either the first element of the list (case 1) or the string itself (case 2)
          if (str(type(scheme))[7:11] == 'list'): # LIST, MANY LEVELS
             outerscheme = scheme[0] # outerscheme equal to the first element in scheme (Case 1)
          else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
             outerscheme = scheme # outerscheme equal to scheme (case 2)


          # THE NUMBER OF FORCEFIELDS PROVIDED MUST EQUAL THE NUMBER
          # OF PROPAGATION LEVELS.  If it does not, we error
          if (len(forcefield) != levels): # if the length of forcefield is not equivalent to levels
             print "[MDL] Error in propagate(): ", levels, " levels of propagation with ", len(forcefield), " force fields."
          outerforcefield = self.deepCopyForce(forcefield[0]) # Take the outermost forcefield and make a deep copy, store it in outerforcefield

          # Now we assemble *chain
          # We have our outermost information (outerscheme, outertime and outerforcefield)
	  # We will call this propagator with the following arguments: outerscheme, outertime, outerforcefield and *chain
	  # In the simplest case (STS), we can simply add params as a fifth parameter and be done
	  # In MTS, we have to add params and THEN add the next integrator in the chain and its params, etc.
          if (str(type(scheme))[7:11] != 'list'): # if the type of scheme is not a list STS (simple case)
             chain += (params,) # chain = chain +(params,) 
          else: # if a list
             if (params.has_key(outerscheme)): # True if params has a key outerscheme
                 chain += (params[outerscheme],)# chain = chain + params[outerscheme]
             else: # if params does not have a key outerscheme - this could potentially be not there (maybe the outer integrator has no parameters)
	           # But, we still must pass an empty list so that the number of parameters is consistent for the propagation
                 chain += ({},) # chain = chain +({},)
 
          # Now for however many levels we have - 1, we add the next propagator in the chain and its parameters
	  for i in range(1, levels): # for i in range from 1, levels
             if (str(type(scheme))[7:11] == 'list' and i < len(scheme)): # if the type of scheme is a list and i is less than the
									 # length of scheme	
                chain += (scheme[i],) # Add the next scheme to the chain
	     # We will need to now add either a cyclelength or a dt (innermost)
	     # We can determine if the integrator is the innermost by comparing i to the length of cyclelength
	     # if i has exceeded this length, we are at the innermost integrator  because that is the only one with no cyclelength
	     # innermost integrator only has dt
             if (str(type(cyclelength))[7:11] == 'list' and i < len(cyclelength)): # if the type of cyclelength is a list and
								            # i is less than the length of cyclelength
                chain += (cyclelength[i],) # Add the cyclelength to the chain
             else: # if not a list
                chain += (dt,) # Add dt to the chain
	     chain += (self.deepCopyForce(forcefield[i]),) # Perform a deep copy of the appropriate force field for this integrator and add to the chain
             if params.has_key(scheme[i]): # Does this integrator have parameters?
                chain += (params[scheme[i]],) # If so, add them to the chain
             else: # if not...
                chain += ({},) # Add an empty list to keep parameter count consistent
       else: #STS, simple.  Four arguments: scheme, dt, forcefield, and parameters (if they are there)
          outertime = dt # set the outertime to dt
          outerscheme = scheme # set outerscheme to scheme
          outerforcefield = self.deepCopyForce(forcefield) # set outerforcefield to deepCopyForce(forcefield)
          chain += (params,) # chain = chain + (params,)
       ############################################################################################################################


       # Build force fields.
       # Tricky, because we could be dealing with
       # a single object or list.
       # In the case of multiple-timestepping (MTS), we evaluate different forces at different timesteps
       # Therefore, the user will provide force fields as a list (i.e., [ff, ff2])
       # In single-timestepping, we would just pass a single forcefield ff for example
 
       # This is the MTS case.  In this case the type of forcefields is a list.  So we have to loop over every element in that list
       if ((str(type(forcefield)))[7:len(str(type(forcefield)))-2] == 'list'): # if the type of forcefield is a list (MTS)
          for ff in forcefield: # step through each element in forcefield
             if (ff.dirty): ff.build() # A dirty bit in a forcefield indicates something has changed, we must rebuild
	     if (ff.gbsa): # Do we want to do Generalized Born?  This is a model for implicit solvent, i.e. simulating a molecule inside a solvent
			   # (i.e. thousands of water molecules) without actually including those molecules.
	        self.phys.myTop.implicitSolvent = 2 # set implicitSolvent to 2 (which means Generalized Born)
		self.phys.myTop.doGBSAOpenMM = 1 # set doGBSAOpenMM to 1 (OpenMM has a generalized Born algorithm, which we use by default)
		self.phys.build() # This requires a rebuilding of the physical object, because we have to change some parameters to account for the
				  # implicit solvent.  The old physical may for instance still assume the molecule is in vacuum.

       # This is the STS case.  Do exactly the above thing, but only on the single forcefield
       else: # if the type of forcefield is not a list 
          if (forcefield.dirty): # if forcefield has a dirty bit
              forcefield.build() # then build forcefield
	  if (forcefield.gbsa): # if forcefield.gbsa
	      self.phys.myTop.implicitSolvent = 2 # set implicitSolvent to 2
              self.phys.myTop.doGBSAOpenMM = 1 # set doGBSAOpenMM to 1
	      self.phys.build() # physical build()


       # Has anything changed in io?  If so, rebuild it.
       # This could happen in the case where we (for example) run a propagator 10 steps, then want to graph something like potential energy
       # Often happens when we want to equilibrate the system for a number of steps before running calculations.
       if (self.io.dirty): # if self.io is using the dirty bit
          self.io.build() # io.build()
        

       
       if (propFactory.getType(outerscheme) == "method"): # if the type fo outerscheme is equal to a "method"
							  # this is true for example in the case of the leapfrog() function
          # Calculate the forces, store them in force.

          ###############################################################
	  # We need a ProtoMolApp, which has been SWIG-wrapped for Python
	  # This portion passes our velocity, positions and forces vectors to that structure
	  # So they then become accessible from C
	  # Note we only need to do this once...
          if (not hasattr(self.phys, 'app')):
             self.phys.app = ProtoMolApp.ProtoMolApp()
             print "Making app 3"
             self.phys.app.makeApp(self.phys.myTop, self.phys.posvec, self.phys.velvec, self.forces.energies, dt)
	  ##############################################################
          self.phys.updateCOM_Momenta()  # Update the center of mass of the system, and atomic momenta (note this equals m*v)
          outerforcefield.calculateForces(self.phys, self.forces)  # Whichever forcefield corresponds to this propagator, calculate forces one time
          self.io.run(self.phys, self.forces, 0, outertime)        # If the user desired I/O, run it.  We may need to plot physical data or forces vs. time
          self.io.myProp = self                                    # The IO structure needs access to me and my data.
          for ii in range(1, steps+1):				   # Run for the appropriate number of steps.
             self.phys.app.energies.clear()			   # Clear energies to zero.
             self.forces.energies.initialize(self.phys)		   # Set energies properly (functions of positions and velocities)
	     # This calls the propagator function with the following arguments:
	     # Physical structure (phys), Forces structure (forces), number of steps (1), current time, corresponding force field.
	     # *chain is empty if we are doing single-timestepping
	     # If we are doing multiple-timestepping, *chain refers to the next propagator function and its arguments
             propFactory.create(1, outerscheme, self.phys, self.forces, self.io, 1, outertime*Constants.invTimeFactor(), outerforcefield, *chain)
	     # Update the physical time with the amount that has passed
             self.phys.time = ii*outertime*Constants.invTimeFactor()
	     # If the user desired I/O, run it.  
	     # Note in the configuration file they set a frequency in steps
	     # In IO.run() we check to see if ii % that frequency is zero, so it will only run that number of steps
             self.io.run(self.phys, self.forces, ii, outertime)
	     # Once the propagator has run, we need to update the center of mass and momenta again.
             self.phys.updateCOM_Momenta()
       else: # Object - this is the case if (1) the integrator is from ProtoMol, or (2) we created a Python class for our integrator that inherited from
	     # ProtoMol's class.  Note that in either case, we are invoking SWIG-wrapped code to some degree as opposed to a pure Python function like
	     # above
	  # This function initializes the propagator
	  # We need the Physical object and the Forces.
          # The create() function returns a new variable to hold the appropriate integrator object.
	  # We pass its scheme name (e.g. "Leapfrog", its outertime which is either cyclelength (MTS) or dt (STS)
	  # Its appropriate force field, and *chain which refers to the next integrator in the chain for MTS or NULL if STS 
	  setPropagator(self, self.phys, self.forces, propFactory.applyModifiers(propFactory.create(1, outerscheme, outertime, outerforcefield, *chain), outerscheme)) # set propagator
          
	  ##################################################################################
	  # This section deals with shake (a bond restraining algorithm)
	  # Shake keeps the bonds near an appropriate specified length
	  shake = False # Shake is turned off by default
          if (params.has_key('shake') and params['shake'] == 'on'): # if params has the key 'shake' and params['shake'] is turned on
              shake = True # We turn shake on
              shakeMod = ModifierShake.ModifierShake(0.000001, 30) # ShakeMod becomes a modifier.  
								   # epsilon is the target amount of deviation from the target length
								   # It is possible that we never hit that amount, so we set maxiterations to 30
              self.myPropagator.adoptPostDriftOrNextModifier(shakeMod) # Shake is run after every update of positions (post-drift)
	  ###################################################################################

	  ###################################################################################
          # This section deals with rattle (also a bond restraining algorithm, but operates on velocities) 
	  rattle = False # Rattle is turned off by default
          if (params.has_key('rattle') and params['rattle'] == 'on'): # if params has key 'rattle' and params['rattle'] is turned on
              rattle = True # We turn rattle on
              rattleMod = self.myPropagator.createRattleModifier(0.02, 30) # Rattle becomes a modifier, epsilon (0.02) is the target deviation
									   # Maximum number of iterations is 30
              self.myPropagator.adoptPostDriftOrNextModifier(rattleMod) # We run Rattle after every update of positions (post-drift)
	  ##################################################################################

          executePropagator(self, self.phys, self.forces, self.io, steps)  # Runs the propagator for a number of steps

	  # Remove the Shake and Rattle modifiers if they were set
	  # This is because we could call prop.propagate a second time, and we do not want to by default run shake and/or rattle again
          if (shake): # if shake
             self.myPropagator.removeModifier(shakeMod) # remove shakemod in myPropagator
          if (rattle): # if rattle
             self.myPropagator.removeModifier(rattleMod)# remove rattlemod in myPropagator

       self.phys.updateCOM_Momenta() # update Center of Mass and angular momentum
           
