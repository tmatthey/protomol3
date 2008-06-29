import Vector3DBlock
import Constants
import numpy
import math
# Leapfrog function accepts initial positions and velocities
# Along with the number of steps to run, the timestep and a group
# of forces
def equilibrate(phys, forces, io, steps, timestep, fg, t0, atomstart):
   """
   Leapfrog propagation method.
   Single timestepping.
   
   @type phys: Physical
   @param phys: The physical system.

   @type forces: Forces
   @param force: MDL Forces object.

   @type io: IO
   @param io: MDL IO object.

   @type steps: int
   @param steps: Number of steps to run.

   @type timestep: float
   @param timestep: Timestep for propagation.

   @type fg: ForceField
   @param fg: MDL force field for evaluation.

   """
   # Update velocities by half a step
   phys.velocities += forces.force*0.5*timestep*phys.invmasses  # half kick
   # Reset solute velocities to zero
   phys.velocities[0:atomstart*3] = 0
   # Update positions by a full step
   phys.positions += phys.velocities*timestep               # drift
   # Calculate new forces with updated position/velocity
   fg.calculateForces(phys, forces)
   # Reset solute forces to zero
   forces.force[0:atomstart*3] = 0
   step = 1
   # Run for the number of passed steps
   while (step < steps):
       # Run I/O
       io.run(phys, forces, step, timestep)
       # Update velocities (full step)
       phys.velocities += forces.force*timestep*phys.invmasses # kick
       # Reset solute velocities to zero
       phys.velocities[0:atomstart*3] = 0
       # Update positions (full step)
       phys.positions += phys.velocities*timestep  # drift
       # Calculate new forces
       fg.calculateForces(phys, forces)
       # Reset solute forces to zero
       forces.force[0:atomstart*3] = 0
       # Scale
       phys.velocities *= math.sqrt(t0/phys.getTemperature())
       # Random
       phys.randomVelocity(t0)
       # Update time
       phys.time = step*timestep       
       # Increment the step
       step = step + 1
   # Update velocities by half a step
   phys.velocities += forces.force*0.5*timestep*phys.invmasses # kick
   # Scale
   phys.randomVelocity(t0)
   # Return positions and velocities as an array of arrays.
   return [phys.positions, phys.velocities]


name="equilibrate"   #: Propagator name for the factory
parameters=('T0', 300, 'atomstart', 0)          #: Parameters and defaults
