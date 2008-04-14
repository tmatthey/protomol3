#include <protomol/integrator/StandardIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/parallel/Parallel.h>

using namespace ProtoMol;
//____ StandardIntegrator
StandardIntegrator::StandardIntegrator() :
  Integrator(), myPreviousIntegrator(NULL) {}

StandardIntegrator::StandardIntegrator(ForceGroup *forceGroup) :
  Integrator(forceGroup), myPreviousIntegrator(NULL) {}

void StandardIntegrator::run(int numTimesteps) {
  for (int i = 0; i < numTimesteps; i++) {
    preStepModify();
    doHalfKick();
    doDriftOrNextIntegrator();
    calculateForces();
    doHalfKick();
    postStepModify();
  }
}

void StandardIntegrator::initializeForces() {
  addModifierBeforeInitialize();
  calculateForces();
  addModifierAfterInitialize();
}

void StandardIntegrator::calculateForces() {
  //  Save current value of potentialEnergy().
  myPotEnergy = app->energies.potentialEnergy();

  myForces->zero();
  preForceModify();

  if (!anyMediForceModify())
    Parallel::distribute(&app->energies, myForces);

  myForcesToEvaluate->evaluateSystemForces(app, myForces);
  mediForceModify();
  myForcesToEvaluate->evaluateExtendedForces(app, myForces);

  if (!anyMediForceModify()) Parallel::reduce(&app->energies, myForces);

  postForceModify();

  //  Compute my potentialEnergy as the difference before/after the call to
  //  calculateForces().
  myPotEnergy = app->energies.potentialEnergy() - myPotEnergy;
}

void StandardIntegrator::doHalfKick() {
  Real h = 0.5 * getTimestep() * Constant::INV_TIMEFACTOR;
  const unsigned int count = app->positions.size();

  updateBeta(h);

  for (unsigned int i = 0; i < count; ++i) {
     app->velocities[i] +=
       (*myForces)[i] * h / app->topology->atoms[i].scaledMass;

  }
  buildMolecularMomentum(&app->velocities, app->topology);
}

void StandardIntegrator::doKick() {
  Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  const unsigned int count = app->positions.size();

  updateBeta(h);

  for (unsigned int i = 0; i < count; ++i)
    app->velocities[i] +=
      (*myForces)[i] * h / app->topology->atoms[i].scaledMass;

  buildMolecularMomentum(&app->velocities, app->topology);
}

Integrator *StandardIntegrator::previous() {
  return myPreviousIntegrator;
}

const Integrator *StandardIntegrator::previous() const {
  return myPreviousIntegrator;
}

