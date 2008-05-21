#include <protomol/integrator/autocorrelator/AutoCorrelatorInner.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;

const string AutoCorrelatorInner::keyword("AutoCorrelatorInner");

AutoCorrelatorInner::AutoCorrelatorInner() :
  LeapfrogIntegrator(), numNonWaters(0) {}

AutoCorrelatorInner::AutoCorrelatorInner(Real timestep,
					 string dcdfile,
					 ForceGroup *overloadedForces) :
  LeapfrogIntegrator(timestep, overloadedForces), numNonWaters(0) {
  myWriter = new DCDTrajectoryWriter(dcdfile, timestep, 0, ISLITTLEENDIAN);
}

AutoCorrelatorInner::~AutoCorrelatorInner() {
  zap(myWriter);
  zap(nonWater);
}

void AutoCorrelatorInner::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  for (unsigned int i = 0; i < myForces->size(); i++)
    if (!app->topology->molecules[app->topology->atoms[i].molecule].water)
      numNonWaters++;
  nonWater = new Vector3DBlock(numNonWaters);
}

void AutoCorrelatorInner::run(int numTimesteps) {

  if (numTimesteps < 1)
    return;
  preStepModify();
  doHalfKickdoDrift();
  calculateForces();
  writeDCD();
  for (int i = 1; i < numTimesteps; i++) {
    doKickdoDrift();
    calculateForces();
    writeDCD();
  }

  doHalfKick();
  postStepModify();

}

void AutoCorrelatorInner::writeDCD()
{
  // Create a structure to hold force
  // values for non-water molecules.
  int pos = 0;
  for (unsigned int i = 0; (i < app->topology->atoms.size()) && pos < numNonWaters; i++) 
    if (!app->topology->molecules[app->topology->atoms[i].molecule].water) {
      (*nonWater)[pos] = (*myForces)[i];
      pos++;
    }   
  // Write the DCD.
  myWriter->write(*nonWater);
}

void AutoCorrelatorInner::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("dcdfile", Value(myDCDFile), ConstraintValueType::NotEmpty()));
}

STSIntegrator *AutoCorrelatorInner::doMake(const vector<Value> &values,
					   ForceGroup *fg) const {
  return new AutoCorrelatorInner(values[0], values[1], fg);
}
