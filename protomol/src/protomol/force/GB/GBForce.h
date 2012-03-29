/* -*- c++ -*- */
#ifndef GBFORCE_H
#define GBFORCE_H

#include <string>

#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  class GBForce {
    public:
      enum {DIST_R2 = 1};
      enum {CUTOFF = 0};
    public:
      GBForce();
      GBForce(Real solute_d, Real solvent_d);
       
      void operator()(Real &energy, Real &force, Real distSquared, Real rDistSquared, 
                      const Vector3D &diff, const GenericTopology *topo, int atom1, 
                      int atom2, ExclusionClass excl) const;
       
      static void accumulateEnergy(ScalarStructure *energies, Real energy);
      static Real getEnergy(const ScalarStructure *energies);
    
      static void postProcess(const GenericTopology *topo, ScalarStructure *energies);
      static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies);

      static bool doParallelPostProcess();
      
      // Parsing
      static std::string getId();
      static unsigned int getParameterSize();
      void getParameters(std::vector<Parameter> &) const;

      static GBForce make(const std::vector<Value> &);
    private:
      Real Force_i_term(const GenericTopology *topo, int atom1) const;
      Real Force_i_j_term(const GenericTopology *topo, const int atom1, const int atom2, const Real ril) const;
    public:
      static const std::string keyword;
    private:
      Real soluteDielec;
      Real solventDielec;
   };
}

#endif
