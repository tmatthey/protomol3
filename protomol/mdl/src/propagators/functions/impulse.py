import Constants
import Vector3DBlock
# Define the impulse integrator function
# It will accept its own parameters for position, velocity,
# number of steps to run, timestep, and force group
# But also accept a next integrator function and its arguments
# because it's multiple time-stepping
def impulse(phys, forces, io, steps, cyclelength, fg, nextinteg, *args):
   """
   Verlet/r-RESPA propagation method.
   Multiple timestepping.
   
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

   @type nextinteg: function
   @param nextinteg: Method handle for next propagator in the chain

   @type args: tuple
   @param args: Parameters for the next propagator in the chain

   """
   # Calculate initial forces
   #fg.calculateForces(phys, forces)
   #fg.calculateForces(phys, phys.positions, phys.velocities, forces.force, forces.energies)
   step = 0
   # For all steps
   timestep = cyclelength*args[0]
   args2 = (args[0]*Constants.invTimeFactor(),)+args[1:len(args)]

   while (step < steps):
      # Update velocities by half a step
      phys.velocities += forces.force*0.5*timestep*phys.invmasses   # half kick
      # Run the next integrator in the chain, and store its results
      # in an array
      nextinteg(phys, forces, io, cyclelength/Constants.invTimeFactor(), *args2)
      # Calculate new forces
      fg.calculateForces(phys, forces)
      # Update velocities by another half step
      phys.velocities += forces.force*0.5*timestep*phys.invmasses   # half kick
      step = step + 1

name="impulse"  #: Propagator name for the factory
parameters=()   #: Parameter names and defaults
     
