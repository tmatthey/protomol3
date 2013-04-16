%module ProtoMolApp
%{
#include <Real.h>
#include <ProtoMolApp.h>
#include <LeapfrogIntegrator.h>
#include <Report.h>
using namespace ProtoMol;
%}

%include <Real.h>
%include <Vector3DBlock.i>
%include <ScalarStructure.i>
%include <ProtoMolApp.h>

%extend ProtoMol::ProtoMolApp {
static void turnOffHints() {
   Report::report << Report::donthint;
}
void uncache() {
   self->outputCache.uncache();
   self->config.registerKeyword("firststep", Value(0));
}
void makeApp(GenericTopology* topo,
             ProtoMol::Vector3DBlock& positions,
             ProtoMol::Vector3DBlock& velocities,
             ScalarStructure energies,
             Real timestep) {
   self->config.registerKeyword("firststep", Value(0));
   self->topology = topo;
   self->positions.vec = positions.vec;
   self->positions.c = positions.c;
   self->velocities.vec = velocities.vec;
   self->velocities.c = velocities.c;
   self->energies = energies;
   self->integrator = new LeapfrogIntegrator(timestep, NULL);
   self->outputCache.initialize(self);
}
};
