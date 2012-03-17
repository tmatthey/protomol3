/* -*- c++ -*- */
#ifndef COULOMBFORCE_H
#define COULOMBFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <string>

namespace ProtoMol {
  //____ CoulombForce
  class CoulombForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom, so the Coulomb
    // constant here is one.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real /*distSquared*/,
                    Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               sqrt(rDistSquared);
      if (excl == EXCLUSION_MODIFIED)
        energy *= topo->coulombScalingFactor;

      force = energy * rDistSquared;
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }
    
    static void postProcess(const GenericTopology *topo, ScalarStructure *energies){
		  
	  }

    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
    }

    //do parallel post process?
    static bool doParallelPostProcess() {
      return false;
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}

    static CoulombForce make(const std::vector<Value> &) {
      return CoulombForce();
    }

  public:
    static const std::string keyword;
  private:
  };
}
#endif /* COULOMBFORCE_H */
