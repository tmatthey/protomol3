#include <protomol/output/Output.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

//____ Output
const string Output::scope("Output");

Output::Output(int freq) :
  app(0), firstStep(0), nextStep(0), lastStep(0), outputFreq(freq) {}

void Output::initialize(const ProtoMolApp *app) {
  this->app = app;

  if (outputFreq <= 0) { // used for finalize only outputs
    if (app->config.valid(InputOutputfreq::keyword)) 
      outputFreq = app->config[InputOutputfreq::keyword];
    else outputFreq = 1;
  }

  if (app->config.valid(InputFirststep::keyword)) {
    nextStep = app->config[InputFirststep::keyword];
    firstStep = app->config[InputFirststep::keyword];
    lastStep = firstStep;
  }

  if (app->config.valid(InputNumsteps::keyword))
    lastStep = lastStep + app->config[InputNumsteps::keyword].operator int();

  doInitialize();
}

void Output::run(int step) {
  if (step >= nextStep) {
    int n = (step - nextStep) / outputFreq;
    nextStep += max(n, 1) * outputFreq;

    if (app->energies.output()) doRun(step);
  }
}

void Output::finalize(int step) {
  if (app->energies.output()) doRun(step);
  doFinalize(step);
}

Output *Output::make(const vector<Value> &values) const {
  assertParameters(values);
  return adjustAlias(doMake(values));
}

bool Output::isIdDefined(const Configuration *config) const {
  if (!addDoKeyword()) return config->valid(getId());

  string str("do" + getId());
  if (!config->valid(getId()) || config->empty(str)) return false;

  if (!config->valid(str)) return true;

  return (bool)(*config)[str];
}
