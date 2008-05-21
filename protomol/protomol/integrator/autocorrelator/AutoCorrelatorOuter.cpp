#include <protomol/integrator/autocorrelator/AutoCorrelatorOuter.h>
#include <protomol/integrator/autocorrelator/AutoCorrelatorInner.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;

const string AutoCorrelatorOuter::keyword("AutoCorrelatorOuter");

AutoCorrelatorOuter::AutoCorrelatorOuter() :
  MTSIntegrator(), numNonWaters(0) {}

AutoCorrelatorOuter::AutoCorrelatorOuter(int cycles,
					 string dcdfile,
					 ForceGroup *overloadedForces, 
					 StandardIntegrator *nextIntegrator) :
  MTSIntegrator(cycles, overloadedForces, nextIntegrator), numNonWaters(0) {
  myWriter = new DCDTrajectoryWriter(dcdfile, getTimestep(), 0, ISLITTLEENDIAN);
}

AutoCorrelatorOuter::~AutoCorrelatorOuter() {
  zap(myWriter);
  zap(nonWater);
}

void AutoCorrelatorOuter::initialize(ProtoMolApp *app) {
  MTSIntegrator::initialize(app);
  // Get number of non-waters here, to save time later
  numNonWaters = ((AutoCorrelatorInner*)myNextIntegrator)->getNumNonWaters();
  nonWater = new Vector3DBlock(numNonWaters);
}

void AutoCorrelatorOuter::run(int numsteps) {
  // Create a structure to hold the current energy values
  ScalarStructure* prevEnergy = new ScalarStructure();
  prevEnergy->intoAssign(app->energies);
  // Calculate forces for this integrator,
  // now the system energies have changed.
  calculateForces();
  // Create a structure to hold force values
  // for non-water molecules.
  int pos = 0;
  for (unsigned int i = 0; (i < app->topology->atoms.size()) && pos < numNonWaters; i++) 
    if (!app->topology->molecules[app->topology->atoms[i].molecule].water) {
      (*nonWater)[pos] = (*myForces)[i];
      pos++;
    }   

  // Write the new forces to the DCD file.
  myWriter->write(*nonWater);
  // Reassign the previous energies, removing
  // the contribution of this force.
  app->energies.intoAssign(*prevEnergy);
  // Clean up, clean up, everybody everywhere...
  delete prevEnergy;
  // Run the next integrator.
  myNextIntegrator->run(numsteps);
}


void AutoCorrelatorOuter::getParameters(vector<Parameter> &parameters)
const {
  MTSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("dcdfile", Value(myDCDFile), ConstraintValueType::NotEmpty()));
}

MTSIntegrator *AutoCorrelatorOuter::doMake(const vector<Value> &values,
					   ForceGroup *fg, 
					   StandardIntegrator *nextIntegrator) const {
  return new AutoCorrelatorOuter(values[0], values[1], fg, nextIntegrator);
}
