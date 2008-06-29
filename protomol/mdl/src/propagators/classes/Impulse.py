from random import *
from MTS import *
from math import *
import Constants

class ImpulseMDL(MTS):

    #  init()  -------------------------------------------------------------  #
    def init(self, phys, forces, prop):
        """
        Initialize propagator: seed the generator, invoke the next
        propagator in the chain, and compute forces.

        @type phys: Physical
        @param phys: The physical system.

        @type forces: Forces
        @param forces: MDL Forces object.

        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        #  Use system time to seed rng.
        prop.initNext(phys, forces)
        prop.calculateForces(forces)

    #  run()  --------------------------------------------------------------  #
    def run(self, phys, forces, prop):
            """
            Run propagator.

            @type phys: Physical
            @param phys: The physical system.

            @type forces: Forces
            @param forces: MDL Forces object.

            @type prop: Propagator
            @param prop: MDL Propagator object.
            """
            phys.velocities += forces.force*0.5*self.cyclelength*prop.myPropagator.next.dt*phys.invmasses   # half kick
            # Run the next integrator in the chain, and store its results
            # in an array
            prop.runNext(phys, forces, self.cyclelength)
            prop.calculateForces(forces)
            phys.velocities += forces.force*0.5*self.cyclelength*prop.myPropagator.next.dt*phys.invmasses

    #  finish()  -----------------------------------------------------------  #
    def finish(self, phys, forces, prop ):
        """
        Finish propagator; in this case just invoke the
        finish method of the next propagator in the chain
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type forces: Forces
        @param forces: MDL Forces object.
        
        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        return
        # Invoked once at simulation finish
        #prop.finishNext(phys, forces, prop)

name="ImpulseMDL"  #: Name of propagation scheme
parameters=()  #: Tuple of parameters
