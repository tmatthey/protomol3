#include <protomol/base/ModuleManager.h>

#include <protomol/module/MainModule.h>
#include <protomol/module/CommandLineModule.h>
#include <protomol/module/ConfigurationModule.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/module/OutputModule.h>
#include <src/module/StringBondedForcesModule.h>
#include <src/module/StringModifierModule.h>
#include <protomol/module/IOModule.h>

#include <protomol/module/IntegratorBaseModule.h>
#include <src/module/StringNormalModeModule.h>
#include <protomol/module/LeapfrogModule.h>
#include <protomol/module/HessianIntegratorModule.h>

#include <protomol/module/NonbondedCutoffForceModule.h>
#include <protomol/module/NonbondedFullForceModule.h>
#include <protomol/module/NonbondedSimpleFullForceModule.h>

using namespace ProtoMol;

void stringModuleInitFunction(ModuleManager *manager) {
  manager->add(new MainModule());
  manager->add(new CommandLineModule());
  manager->add(new ConfigurationModule());
  manager->add(new TopologyModule());
  manager->add(new OutputModule());
  manager->add(new StringModifierModule());
  manager->add(new IOModule());

  // Integrators
  manager->add(new IntegratorBaseModule());
  manager->add(new StringNormalModeModule());
  manager->add(new LeapfrogModule());
  manager->add(new HessianIntegratorModule());

  // Forces
  manager->add(new StringBondedForcesModule());
  manager->add(new NonbondedCutoffForceModule());
  manager->add(new NonbondedFullForceModule());
  manager->add(new NonbondedSimpleFullForceModule());
}

