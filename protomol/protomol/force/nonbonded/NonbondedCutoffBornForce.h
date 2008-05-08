/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFBORNFORCE_H
#define NONBONDEDCUTOFFBORNFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffForce.h>
#include <protomol/force/nonbonded/Born.h>
#include <iostream>
#include <iomanip>

namespace ProtoMol {
  //____ NonbondedCutoffSystemForce

  template<class TCellManager, class TOneAtomPair>
  class NonbondedCutoffBornForce :
    public NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                                NonbondedCutoffBornForce<TCellManager,
                                                         TOneAtomPair> >,
    public Born<TOneAtomPair> {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    typedef typename TOneAtomPair::BoundaryConditions BoundaryConditions;
    typedef Topology<BoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef typename ProtoMol::CellPair CellPairType;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedCutoffBornForce() {}
    NonbondedCutoffBornForce(Real cutoff, TOneAtomPair oneAtomPair) :
      NonbondedCutoffForce<TCellManager, TOneAtomPair, SystemForce,
                           NonbondedCutoffBornForce>(cutoff, oneAtomPair) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(GenericTopology *topo, const Vector3DBlock *positions,
                          Vector3DBlock *forces, ScalarStructure *energies) {
      // Standard for forces
      if (!topo->doSCPISM)
        THROW("To use SCPISM forces, please set doscpism in the configuration "
              "file");
      
      RealTopologyType *realTopo = dynamic_cast<RealTopologyType *>(topo);
      this->myOneAtomPair.getNonbondedForceFunction()->
        resize(topo->atoms.size());
      for (unsigned int i = 0; i < topo->atoms.size(); i++)
        this->myOneAtomPair.getNonbondedForceFunction()->
          resizeDR(i, topo->atoms.size());
      
      this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
      realTopo->updateCellLists(positions);
      this->enumerator.initialize(realTopo, this->myCutoff);
      this->doEvaluate(topo, realTopo->cellLists.size());
      evaluateBorn(this->myOneAtomPair, topo, forces, energies);
    }

    virtual void parallelEvaluate(GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies) {
      RealTopologyType *realTopo = dynamic_cast<RealTopologyType *>(topo);
      
      this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
      realTopo->updateCellLists(positions);
      this->enumerator.initialize(realTopo, this->myCutoff);
      
      unsigned int n = realTopo->cellLists.size();
      unsigned int count = numberOfBlocks(realTopo, positions);
      
      for (unsigned int i = 0; i < count; i++) {
        unsigned int l = (n * (i + 1)) / count - (n * i) / count;

        if (Parallel::next()) this->doEvaluate(topo, l);
        else this->enumerator.nextNewPair(l);
      }
    }

    virtual std::string getKeyword() const {return "NonbondedCutoffBorn";}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Vector3DBlock bforces;
  };
}
#endif /* NONBONDEDCUTOFFSYSTEMFORCE_H */
