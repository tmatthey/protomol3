from random import *
from STS import *
from math import *
import Constants

class LeapfrogMDL(STS):
    def init(self, phys, forces, prop):
        # Update velocities by half a step
        phys.velocities += forces.force*0.5*self.dt*phys.invmasses  # half kick
        # Update positions by a full step
        phys.positions += phys.velocities*self.dt               # drift
        # Calculate new forces with updated position/velocity
        prop.calculateForces(forces)


    def run(self, phys, forces, prop):
        # Update velocities (full step)
        phys.velocities += forces.force*self.dt*phys.invmasses # kick
        # Update positions (full step)
        phys.positions += phys.velocities*self.dt  # drift
        # Calculate new forces
        prop.calculateForces(forces)

    def finish(self, phys, forces, prop):
        # Update velocities by half a step
        phys.velocities += forces.force*0.5*self.dt*phys.invmasses # kick

name="LeapfrogMDL"  #: Name of propagation scheme
parameters=() #: Parameter names and defaults
