import Vector3DBlock
import Constants
import numpy
# Simplified Takahashi function accepts initial positions and velocities
# Along with the number of steps to run, the timestep and a group
# of forces
def takahashi(phys, forces, io, steps, timestep, fg):
   """
   Simplified Takahashi-Imada propagation method.
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
   # Update positions by a full step
   phys.positions += phys.velocities*timestep               # drift
   # Calculate new forces with updated position/velocity
   fg.calculateForces(phys, forces)
   # update force argument
   auxpositions = phys.positions.copy()
   phys.positions +=  timestep * timestep / 12.0 * phys.invmasses * forces.force
   fg.calculateForces(phys, forces)

   phys.positions.put(range(0, 3*phys.N()), auxpositions)


   step = 1
   # Run for the number of passed steps
   while (step < steps):
       # Run I/O
       io.run(phys, forces, step, timestep)
       # Update velocities (full step)
       phys.velocities += forces.force*timestep*phys.invmasses # kick
       # Update positions (full step)
       phys.positions += phys.velocities*timestep  # drift
       # Calculate new forces
       fg.calculateForces(phys, forces)
       auxpositions.put(range(0, 3*phys.N()), phys.positions)
       phys.positions +=  timestep * timestep / 12.0 * phys.invmasses * forces.force
       fg.calculateForces(phys, forces)

       phys.positions.put(range(0, 3*phys.N()), auxpositions)

       # Update time
       phys.time = step*timestep       
       # Increment the step
       step = step + 1
   # Update velocities by half a step
   phys.velocities += forces.force*0.5*timestep*phys.invmasses # kick
   # Return positions and velocities as an array of arrays.
   return [phys.positions, phys.velocities]


name="takahashi"   #: Propagator name for the factory
parameters=()      #: Parameter names and defaults
