import numpy
import math
#from Potential import *
from PySystemForce import *
import sys
import Constants
def norm(v3d):
    return numpy.sqrt(v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2])

def norm2(v3d):
    return v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2]
  
class ElectrostaticForce:
  """
  Implement a harmonic dihedral constraining potential:
  U(x) = k*(phi - phi0)^2
  """
  
  def __init__(self, phys, forces):
    """
    Initialize an object of type HDForce
    
    @type phys: Physical
    @param phys: The physical system.

    @type forces: Forces
    @param forces: MDL Forces object

    """
    self.phys = phys
    self.forces = forces

    n = self.phys.numAtoms()
    self.q = []
    
    for i in range (0, n):
        self.q.append(self.phys.charge(i+1))
      
  def eval(self):
    """
    Modify energy and force vector to include this force term.
    """
    n = self.phys.numAtoms()
    print "CE Before: ", self.forces.energies.coulombEnergy(self.phys)
#    values = []
    for i in range (0, n):
      for j in range (i+1, n):
#        if (values.count("(%(i)d,%(j)d)" % {'i': i, 'j': j}) > 0):
#          print "(%(i)d,%(j)d)" % {'i': i, 'j': j}
#        else:
#          values.append("(%(i)d,%(j)d)" % {'i': i, 'j': j})
        rij = self.phys.positions[j*3:j*3+3] - self.phys.positions[i*3:i*3+3]
        rdistsquared = 1.0 / norm2(rij)
        same = (self.phys.myTop.getAtom(i).molecule == self.phys.myTop.getAtom(j).molecule)
        en = self.q[i]*self.q[j]*numpy.sqrt(rdistsquared)*self.phys.myTop.checkExclusion(i,j,same)
#        if (i == 5 and j == 253):
#          print "(", i, ",", j, ")"
#          print "  Coords for ", i, ": ", self.phys.positions[i*3:i*3+3]
#          print "  Coords for ", j, ": ", self.phys.positions[j*3:j*3+3]
#          print "  Energy: ", en
#          print "  Scaled Charge i: ", self.q[i]
#          print "  Scaled Charge j: ", self.q[j]
#          print "  distSquared: ", norm2(rij)
#          print "  rDistSquared: ", rdistsquared
#        if (self.phys.myTop.checkExclusion(i,j,same) != 0):
#          file = open('MDLab Pairs.txt', 'a')
#          file.write("(%(i)d,%(j)d) : %(en)f\n" % {'i': i, 'j': j, 'en': en})
#          file.close()
#          print "Energy for atom pair ", i, ", ", j, ": ", en, " | Exclusion: ", self.phys.myTop.checkExclusion(i,j,same)
#          print "Atom: ", i, ", ", j
#          print "Python Energy: ", en
#          print "Python rDistSquared: ", rdistsquared
#          print "Python SQRT rDistSquared: ", numpy.sqrt(rdistsquared)
#          print "Python Charges: ", self.q[i], ", ", self.q[j]
#          print "Excluded: ", self.phys.myTop.checkExclusion(i,j)
        self.forces.energies.addCoulombEnergy(en, self.phys)
        fo = en*rdistsquared
#        if (self.phys.myTop.checkExclusion(i,j,same) != 0):
#          file = open('MDLab Force Values.txt', 'a')
#          file.write("(%(i)d,%(j)d) : %(fo)f\n" % {'i': i, 'j': j, 'fo': fo})
#          file.close()
        self.forces.force[i*3:i*3+3] -= fo * rij
        self.forces.force[j*3:j*3+3] += fo * rij
    print "CE After: ", self.forces.energies.coulombEnergy(self.phys)
