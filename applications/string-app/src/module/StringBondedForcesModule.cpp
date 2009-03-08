#include <src/module/StringBondedForcesModule.h>

#include <protomol/force/bonded/DihedralSystemForce.h>
#include <protomol/force/bonded/BondSystemForce.h>
#include <protomol/force/bonded/AngleSystemForce.h>
#include <protomol/force/bonded/ImproperSystemForce.h>
#include <protomol/force/bonded/HarmDihedralSystemForce.h>
#include <src/force/bonded/HarmNMRestSystemForce.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

using namespace std;
using namespace ProtoMol;

void StringBondedForcesModule::init(ProtoMolApp *app) {}

void StringBondedForcesModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;

  string boundConds =
    app->config[InputBoundaryConditions::keyword];

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    f.registerExemplar(new DihedralSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new BondSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new AngleSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new ImproperSystemForce<PeriodicBoundaryConditions>());
    f.registerExemplar(new
                       HarmDihedralSystemForce<PeriodicBoundaryConditions>());

  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    f.registerExemplar(new DihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new BondSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new AngleSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new ImproperSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new HarmDihedralSystemForce<VacuumBoundaryConditions>());
    f.registerExemplar(new HarmNMRestSystemForce<VacuumBoundaryConditions>());
  }
}
