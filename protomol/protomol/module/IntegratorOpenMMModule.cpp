#include <protomol/module/IntegratorOpenMMModule.h>

#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/integrator/openMM/NormalModeOpenMM.h>

#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

void IntegratorOpenMMModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new OpenMMIntegrator());
  app->integratorFactory.registerExemplar(new NormalModeOpenMM());

}
