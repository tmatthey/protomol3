import TopologyUtilities
from ScalarStructure import *

######################################################################
# Class Energies
# Function: 
#####################################################################

class Energies(ScalarStructure):
   """
   Holds the system energies
   These include: potential (bond, angle, ...) and kinetic.

   Inherits from precompiled class ScalarStructure
   This manages energies in a table.  The inherited method getTable()
   will index this table.
   
   This provides several inherited methods which are used:
     - potentialEnergy()
     - intoAdd()
     - intoSubtract()
     - molecularVirial()
     - computeVirial()
   
   """
   def initialize(self, phys):
      self.phys = phys #: Physical object
      
   def computeMolecularVirial(self):
     # the Virial theorem is used to describe the average Kinetic Energy over a period of time. 
     # consisting of a stable system of N particles or molecules.
     # Tensors are ways of describing a system of vectors or scalars in relation to objects (like molecules in this case)
     # Tensors include dot products -- Can be expressed in two different ways: Algebraically and Geometrically
     # Algebraically it is the sum of the products of the corresponding vectors.
     # Geometrically it is the magnitude of the corresponding vectors and the cosine of the angle between them resulting in a scalar projection
     # Cross product is a way of calculating the direction of the vector perpendicular to the original vectors
     # Cross product Geometrically is the magnitude of A and B multiplied by sine of the angle.
     # this results in vectors being either up, down, or zero depending on the direction of the vector and if they are perpendicular
     # if they are parallel it will result in a zero vector.
      """
      Tells the energies structure to include the molecular
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's molecularVirial()
      return self.phys.app.energies.molecularVirial(1)

   def computeVirial(self):
      """
      Tells the energies structure to include the 
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's virial()
      return self.phys.app.energies.virial(1)


   def addBondEnergy(self, r):
      """
      Accumulate into the bond energy
      
      @type r: float
      @param r: Quantity to accumulate.
      
      """
      self.phys.app.energies.setTable(2, self.bondEnergy()+r) # this sets up a table to acquire bond energy values into a table
      # at position 2

   def addAngleEnergy(self, r):
      """            
      Accumulate into the angle energy

      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(3, self.angleEnergy()+r) # this sets up a table to acquire Angle energy values into a table
      # at position 3      

   def addDihedralEnergy(self, r):
      """
      Accumulate into the dihedral energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(4, self.dihedralEnergy()+r) # this sets a table to acquire Dihedral energy values into a table
      # at position 4

   def addImproperEnergy(self, r):
      """
      Accumulate into the improper energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(5, self.improperEnergy()+r) # this sets a table to acquire Improper energy values into a table
      # at position 5

   def addShadowEnergy(self, r):
      """
      Accumulate into the shadow energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(34, self.shadowEnergy()+r) # this sets a table to acquire Shadow energy values into a table
      # at position 34

   def addCoulombEnergy(self, r):
      """
      Accumulate into the electrostatic energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(0, self.coulombEnergy()+r) # this sets a table to acquire Coulomb energy values into a table 
      # at position 0

   def addLJEnergy(self, r):
      """
      Accumulate into the van der Waals energy

      @type r: float
      @param r: Quantity to accumulate.      
      """
      self.phys.app.energies.setTable(1, self.ljEnergy()+r) # this sets a table to acquire Lennard Jones energy values into a table
      # at position 1

            
   def coulombEnergy(self, phys):
      """
      @rtype: float
      @return: Electrostatic energy
      """

      # Coulomb energy is table index 1
      return phys.app.energies.getTable(0) # returns what is at positon 0 in the table, which is Coulomb Energy
   
   def ljEnergy(self, phys):
      """
      @rtype: float
      @return: van der Waals energy
      """

      # LJ energy is table index 2
      return phys.app.energies.getTable(1) # returns what is at positon 1 in the table, which is Lennard Jones energy

   def bondEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to two-atom bond deviations from equilibrium lengths.
      """

      # Bond energy is table index 3
      return phys.app.energies.getTable(2) # returns what is at positon 2 in the table, which is Bond Energy

   def angleEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to three-atom angle deviations from equilibrium angles.
      """

      # Angle energy is table index 4
      return phys.app.energies.getTable(3)# returns what is at positon 3 in the table, which is Angle Energy

   def dihedralEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to four-atom dihedral deviations from equilibrium angles.
      """

      # Dihedral energy is table index 5
      return phys.app.energies.getTable(4) # returns what is at positon 4 in the table, which is Dihedral Energy
   
   def improperEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to four-atom improper deviations from equilibrium angles.
      """

      # Improper energy is table index 6
      return phys.app.energies.getTable(5) # returns what is at positon 5 in the table, which is Improper Energy

   def shadowEnergy(self, phys):
      """
      @rtype: float
      @return: Shadow energy.
      """

      # Shadow energy is table index 4
      return phys.app.energies.getTable(34) # returns what is at positon 34 in the table, which is shadow Energy

   def kineticEnergy(self, phys):
    # Kinetic energy: .5mv^2
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Kinetic energy, as a sum of 0.5*m*v^2 for each atom.
      """
      return TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)

   def potentialEnergy(self, phys):
    # potential energy is equivalent to mass x gravity x height
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Potential energy
      """
      return phys.app.energies.potentialEnergy()

   def totalEnergy(self, phys):
    # total energy is equal to the addition of Kinetic and Potential
    # Kinetic energy: .5*m*v^2
    # Potential energy: mgh
    # therefore total energy: .5mv^2 + mgh
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Total energy, as a sum of potential and kinetic.
      """

      # potentialEnergy() is a member function of ScalarStructure
      return phys.app.energies.potentialEnergy()+TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)
   
