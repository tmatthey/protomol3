%module STSIntegrator
%{
#include "STSIntegrator.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/integrator/StandardIntegrator.h>
#include <protomol/integrator/STSIntegrator.h>
#include "ndarrayobject.h"
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%feature ("dynamic_cast");
%feature("notabstract") ProtoMol::STSIntegrator;
%include "std_string.i"
%include <protomol/integrator/Integrator.h>
%include <protomol/integrator/StandardIntegrator.h>
%include <protomol/integrator/STSIntegrator.h>
//%include "STSIntegrator.h"


%extend ProtoMol::STSIntegrator {
ProtoMolApp* appInit(GenericTopology* topo,
                     Vector3DBlock& positions,
                      Vector3DBlock& velocities,
                      ScalarStructure energies) {
   import_array1(NULL);
   ProtoMolApp* app = new ProtoMolApp();
   app->topology = topo;
   app->positions.vec = positions.vec;
   app->positions.c = positions.c;
   app->velocities.vec = velocities.vec;
   app->velocities.c = velocities.c;
   app->energies = energies;
   self->initialize(app);
   app->integrator = self;
   app->outputCache.initialize(app);
   return app;
}
};
