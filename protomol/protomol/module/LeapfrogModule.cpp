#include <protomol/module/LeapfrogModule.h>

#include <protomol/ProtoMolApp.h>

#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/LeapfrogTruncatedShadow.h>
#include <protomol/integrator/leapfrog/DMDLeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/PLeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/NoseNVTLeapfrogIntegrator.h>

// Autocorrelators are added here because they run Leapfrog internally
#include <protomol/integrator/autocorrelator/AutoCorrelatorOuter.h>
#include <protomol/integrator/autocorrelator/AutoCorrelatorInner.h>

using namespace std;
using namespace ProtoMol;

void LeapfrogModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new LeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new LeapfrogTruncatedShadow());
  app->integratorFactory.registerExemplar(new DMDLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new PLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new NoseNVTLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new AutoCorrelatorOuter());
  app->integratorFactory.registerExemplar(new AutoCorrelatorInner());
}
