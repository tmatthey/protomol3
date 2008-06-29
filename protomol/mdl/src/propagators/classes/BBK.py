from random import *
from STS import *
from math import *
import Constants
import _TopologyUtilities
import numpy

class BBK(STS):
    def init(self, phys, forces, prop):
       self.gamma = self.gamma*0.001/Constants.invTimeFactor()
       prop.calculateForces(forces)

    def run(self, phys, forces, prop):
       forceconstant = 2*Constants.boltzmann()*self.temp*self.gamma/self.dt      # assign random force
       forces.force += forces.randomForce(phys,self.seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities *= (1.0-0.5*self.dt*self.gamma)         # first half kick
       phys.velocities += forces.force*0.5*self.dt*phys.invmasses
       #_topologyutilities.buildMolecularMomentum(phys.velvec,phys.myTop)
       phys.positions += phys.velocities*self.dt                # drift
       #_topologyutilities.buildMolecularCenterOfMass(phys.posvec,phys.myTop)
       prop.calculateForces(forces)
       forceconstant = 2*Constants.boltzmann()*self.temp*self.gamma/self.dt      # assign random force
       forces.force += forces.randomForce(phys,self.seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities += forces.force*0.5*self.dt*phys.invmasses   # second first kick
       phys.velocities *= (1.0/(1.0+0.5*self.dt*self.gamma))
       #_topologyutilities.buildMolecularMomentum(phys.velvec,phys.myTop)
       #prop.calculateForces(forces)

name="BBK"  #: Name of propagation scheme
parameters=("temp", 300,
            "gamma", 2,
            "seed", 1234) #: Parameter names and defaults
