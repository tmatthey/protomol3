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
   def __add__(self, rhs):
      """
      Overloads the '+' operator to add two energy structures
      This simply adds corresponding energy terms for cumulative
      bond energy, etc.

      @type rhs: Energies
      @param rhs: The second Energies object.

      @rtype: Energies
      @return: The sum of the two Energies objects.
      """

      # Invoke ScalarStructure's intoAdd()
      retval = self
      retval.intoAdd(rhs)
      return retval


   def __sub__(self, rhs):
      """
      Overloads the '-' operator to add two energy structures
      This simply subtracts corresponding energy terms for cumulative
      bond energy, etc.

      @type rhs: Energies
      @param rhs: The second Energies object.

      @rtype: Energies
      @return: The difference of the two Energies objects.
      """
      retval = self
      # Invoke ScalarStructure's intoSubtract()
      retval.intoSubtract(rhs)
      return retval


   def computeMolecularVirial(self):
      """
      Tells the energies structure to include the molecular
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's molecularVirial()
      return self.molecularVirial(1)


   def computeVirial(self):
      """
      Tells the energies structure to include the 
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's virial()
      return self.virial(1)


   def addBondEnergy(self, r):
      """
      Accumulate into the bond energy
      
      @type r: float
      @param r: Quantity to accumulate.
      
      """
      self.setTable(3, self.bondEnergy()+r)

   def addAngleEnergy(self, r):
      """            
      Accumulate into the angle energy

      @type r: float
      @param r: Quantity to accumulate.
      """
      self.setTable(4, self.angleEnergy()+r)
      
   def addDihedralEnergy(self, r):
      """
      Accumulate into the dihedral energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.setTable(5, self.dihedralEnergy()+r)
      #self.setTable(5, self.dihedralEnergy()+r)

   def addImproperEnergy(self, r):
      """
      Accumulate into the improper energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.setTable(6, self.improperEnergy()+r)

   def addShadowEnergy(self, r):
      """
      Accumulate into the shadow energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.setTable(34, self.shadowEnergy()+r)

   def addCoulombEnergy(self, r):
      """
      Accumulate into the electrostatic energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.setTable(1, self.coulombEnergy()+r)
      #self.setTable(1, 4)

   def addLJEnergy(self, r):
      """
      Accumulate into the van der Waals energy

      @type r: float
      @param r: Quantity to accumulate.      
      """
      self.setTable(2, self.ljEnergy()+r)

            
   def coulombEnergy(self):
      """
      @rtype: float
      @return: Electrostatic energy
      """

      # Coulomb energy is table index 1
      return self.getTable(1)
   
   def ljEnergy(self):
      """
      @rtype: float
      @return: van der Waals energy
      """

      # LJ energy is table index 2
      return self.getTable(2)

   def bondEnergy(self):
      """
      @rtype: float
      @return: Energy due to two-atom bond deviations from equilibrium lengths.
      """

      # Bond energy is table index 3
      return self.getTable(3)

   def angleEnergy(self):
      """
      @rtype: float
      @return: Energy due to three-atom angle deviations from equilibrium angles.
      """

      # Angle energy is table index 4
      return self.getTable(4)

   def dihedralEnergy(self):
      """
      @rtype: float
      @return: Energy due to four-atom dihedral deviations from equilibrium angles.
      """

      # Dihedral energy is table index 5
      return self.getTable(5)
   
   def improperEnergy(self):
      """
      @rtype: float
      @return: Energy due to four-atom improper deviations from equilibrium angles.
      """

      # Improper energy is table index 6
      return self.getTable(6)

   def shadowEnergy(self):
      """
      @rtype: float
      @return: Shadow energy.
      """

      # Shadow energy is table index 4
      return self.getTable(34)

   def kineticEnergy(self, phys):
      """
      @rtype: float
      @return: Kinetic energy, as a sum of 0.5*m*v^2 for each atom.
      """
      return TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)

   def totalEnergy(self, phys):
      """
      @rtype: float
      @return: Total energy, as a sum of potential and kinetic.
      """

      # potentialEnergy() is a member function of ScalarStructure
      return self.potentialEnergy()+TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)
   
