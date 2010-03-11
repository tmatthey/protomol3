#include "VesselForcesModule.h"

#include <protomol/force/bonded/RBDihedralSystemForce.h>
#include <protomol/force/bonded/DihedralSystemForce.h>
#include "../force/VesselForce.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

using namespace std;
using namespace ProtoMol;

void VesselForcesModule::init(ProtoMolApp *app) {}

void VesselForcesModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;

  string boundConds =
    app->config[InputBoundaryConditions::keyword];

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
    f.registerExemplar(new VesselForce<PeriodicBoundaryConditions>());

  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
    f.registerExemplar(new VesselForce<VacuumBoundaryConditions>());

  }
}
