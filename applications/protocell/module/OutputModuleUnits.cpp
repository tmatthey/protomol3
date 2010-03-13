#include "OutputModuleUnits.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/config/Configuration.h>
#include <protomol/factory/OutputFactory.h>

#include "../output/OutputScreenUnit.h"

using namespace std;
using namespace ProtoMol;

void OutputModuleUnits::init(ProtoMolApp *app) {
  OutputFactory &f = app->outputFactory;

  f.unregisterExemplar("Screen");
  
  f.registerExemplar(new OutputScreenUnit());

}
