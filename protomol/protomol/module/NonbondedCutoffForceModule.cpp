#include <protomol/module/NonbondedCutoffForceModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

#include <protomol/force/OneAtomPair.h>
#include <protomol/force/OneAtomPairTwo.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/coulomb/CoulombBornRadiiForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffSystemForce.h>
#include <protomol/force/nonbonded/NonbondedCutoffBornForce.h>
#include <protomol/force/table/LennardJonesTableForce.h>

#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/UniversalSwitchingFunction.h>

#include <protomol/topology/CellListEnumeratorPeriodicBoundaries.h>
#include <protomol/topology/CellListEnumeratorStandard.h>

using namespace std;
using namespace ProtoMol;

void NonbondedCutoffForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  // To make this a bit more readable
  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef CubicCellManager CCM;
  typedef C1SwitchingFunction C1;
  typedef C2SwitchingFunction C2;
  typedef CnSwitchingFunction Cn;
  typedef CmpCnCnSwitchingFunction CmpCnCn;
  typedef UniversalSwitchingFunction Universal;
#define CutoffSystem NonbondedCutoffSystemForce
#define CutoffBorn NonbondedCutoffBornForce

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    // NonbondedCutoffSystemForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C1, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C2, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, CmpCnCn, CoulombForce> >());
    
    // NonbondedCutoffSystemForce LennardJonesForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C1, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, C2, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Cn, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, CmpCnCn,
          LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<C2, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<Cn, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<PBC, Universal,
          LennardJonesTableForce<CmpCnCn, 7 ,Real> > >());
    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce,
          C1, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, C2, LennardJonesForce,
          Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, Cn, LennardJonesForce,
          Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<PBC, CmpCnCn, LennardJonesForce,
          CmpCnCn, CoulombForce> >());


  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    // NonbondedCutoffSystemForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn, CoulombForce> >());

    // NonbondedCutoffSystemForce LennardJonesForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn,
          LennardJonesForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<C2, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<Cn, 7 ,Real> > >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Universal,
          LennardJonesTableForce<CmpCnCn, 7 ,Real> > >());
    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C1,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, Cn,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, Cn, LennardJonesForce, Cn,
          CoulombForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, CmpCnCn, LennardJonesForce,
          C1, CoulombForce> >());

    // SCPISM stuff
    // NonbondedCutoffSystemForce CoulombSCPISMForce
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C1, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, C2, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, Cn, CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPair<VBC, CmpCnCn,
          CoulombSCPISMForce> >());

    // OneAtomPairTwo
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C1,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, C2, LennardJonesForce, C2,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, Cn, LennardJonesForce, Cn,
          CoulombSCPISMForce> >());
    f.reg(new CutoffSystem<CCM, OneAtomPairTwo<VBC, CmpCnCn, LennardJonesForce,
          C1, CoulombSCPISMForce> >());

    // CutoffBorn CoulombSCPISMForce
    f.reg(new CutoffBorn<CCM, OneAtomPair<VBC, C1, CoulombBornRadiiForce> >());
    f.reg(new CutoffBorn<CCM, OneAtomPair<VBC, C2, CoulombBornRadiiForce> >());
    f.reg(new CutoffBorn<CCM, OneAtomPair<VBC, Cn, CoulombBornRadiiForce> >());
    f.reg(new CutoffBorn<CCM, OneAtomPair<VBC, CmpCnCn,
          CoulombBornRadiiForce> >());
  }
}
