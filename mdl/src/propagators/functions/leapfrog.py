import Vector3DBlock
import Constants
import numpy
# Leapfrog function accepts initial positions and velocities
# Along with the number of steps to run, the timestep and a group
# of forces
def leapfrog(phys, forces, io, steps, timestep, fg):
   """
   Implements the Leapfrog method.
      1. Half-timestep update of velocities.
      2. Full-timestep update of positions.
      3. Half-timestep update of velocities.
   cf. R. W. Hockney and J. W. Eastwood, Computer Simulation Using Particles.
   New York: McGraw-Hill, 1981.
   
   @type phys: Physical
   @param phys: The physical system.

   @type forces: Forces
   @param forces: MDL Forces object.

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
#   print (phys.invmasses*forces.force*0.5).matrix
#   print "Invmasses Before:"
#   print phys.gpuInvmasses.get_matrix()
#   raw_input()
#   print "TEST"
#   print (phys.invmasses*0.5*timestep).get_matrix()#*0.5*timestep).get_matrix()
#   raw_input()
#   print "Forces: "
#   print forces.force.get_matrix()
#   raw_input()
#   print "Value to add:"
#   print (phys.invmasses*forces.force*0.5*timestep).get_matrix()
#   raw_input()
#   print "Velocities before change:"
#   print phys.gpuVelocities.get_matrix()
#   raw_input()
#   print "Forces:"
#   print forces.forcevec.getC()
#   raw_input()
#   print "Result:"
#   print (phys.invmasses*forces.force*0.5*timestep).get_matrix()
#   raw_input()
#   phys.invmasses*forces.force*0.5*timestep
#   print "Invmasses After:"
#   print phys.gpuInvmasses.get_matrix()
#   raw_input()
   phys.velocities += phys.invmasses*forces.force*0.5*timestep  # half kick
#   print "Invmasses After:"
#   print phys.gpuInvmasses.get_matrix()
#   raw_input()
#   print "Velocities after change:"
#   print phys.gpuVelocities.get_matrix()
#   raw_input()
#   print "0: ", phys.positions.get_matrix()
#   raw_input()
#   print "VEL 1: ", phys.velocities.matrix
#   raw_input()
#   print "CVEL 1: ", phys.velvec.getC()
   # Update positions by a full step
#   print "CPOS: ", phys.posvec.getC()
#   raw_input()
#   print (phys.positions + (phys.velocities * timestep)).get_matrix()
#   raw_input()
#   print "Vel: ", phys.velocities.get_matrix()
#   raw_input()
#   print "Pos before change:"
#   print phys.gpuPositions.get_matrix()
#   raw_input()
   phys.positions += phys.velocities*timestep               # drift 
#   print "Pos after change:"
#   print phys.gpuPositions.get_matrix()
#   raw_input()
#   print "1: ", phys.positions.get_matrix()
#   raw_input()
#   print "POSITIONS: ", phys.positions.get_matrix()
#   raw_input()
#   print "Pos: ", phys.positions[0]
#   print "CPOS: ", phys.posvec.getC()
#   raw_input()
#   print "Forces before change:"
#   print forces.gpuForce.get_matrix()
#   raw_input()
   # Calculate new forces with updated position/velocity
   fg.calculateForces(phys, forces)
#   print "Forces after change:"
#   print forces.gpuForce.get_matrix()
#   raw_input()
#   if (phys.gpu):
#     print "Forces: "
#     print forces.force.get_matrix()
#   raw_input()
#   print "2: ", phys.positions.get_matrix()
#   raw_input()
   step = 1
   # Run for the number of passed steps
   while (step < steps):
       # Run I/O
       io.run(phys, forces, step, timestep)
       # Update velocities (full step)
       phys.velocities += phys.invmasses*forces.force*timestep#*phys.invmasses # kick
       # Update positions (full step)
       phys.positions += phys.velocities*timestep  # drift
       # Calculate new forces
       fg.calculateForces(phys, forces)
       # Update time
       phys.time = step*timestep       
       # Increment the step
       step = step + 1
   # Update velocities by half a step
#   print "VEL 2: ", phys.velocities[0]
#   print "CVEL 2: ", phys.velvec.getC()[0]
#   print "INVMASSES: ", phys.invmasses[0]
#   print "FORCES: ", forces.force[0]
#   print "F: ", forces.force[0]
#   print (phys.invmasses).get_matrix()
#   print "Forces at Second Change:"
#   print forces.gpuForce.get_matrix()
#   raw_input()
#   print "Velocities Before Second Change:"
#   print phys.gpuVelocities.get_matrix()
#   raw_input()
   phys.velocities += phys.invmasses*forces.force*0.5*timestep#*phys.invmasses # half kick
#   print "Velocities After Second Change:"
#   print phys.gpuVelocities.get_matrix()
#   raw_input()
#   print "3: ", phys.positions.get_matrix()
#   raw_input()
#   print "VEL 3: ", phys.velocities[0]
#   print "CVEL 3: ", phys.velvec.getC()[0]
#   print phys.positions.matrix
#   phys.velvec.setC(phys.velocities.matrix)
#   phys.posvec.setC(phys.positions.matrix)
   # Return positions and velocities as an array of arrays.
   return [phys.positions, phys.velocities]


name="leapfrog"   #: Propagator name for the factory
parameters=()     #: Parameters and defaults
