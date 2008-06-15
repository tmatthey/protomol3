%module ProtoMolApp
%{
#include <protomol/ProtoMolApp.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/type/Vector3DBlock.i>
%include <protomol/type/ScalarStructure.i>
%include <protomol/ProtoMolApp.h>

%extend ProtoMol::ProtoMolApp {
void uncache() {
   self->outputCache.uncache();
}
void makeApp(GenericTopology* topo,
             ProtoMol::Vector3DBlock& positions,
             ProtoMol::Vector3DBlock& velocities,
             ScalarStructure energies,
             Real timestep) {
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
