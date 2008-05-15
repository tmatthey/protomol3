%module ProtoMolApp
%{
#include "ProtoMolApp.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/STSIntegrator.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/type/Vector3DBlock.i>
%include <protomol/type/ScalarStructure.h>
%include "ProtoMolApp.h"

%extend ProtoMol::ProtoMolApp {
void makeApp(GenericTopology* topo,
             ProtoMol::Vector3DBlock positions,
             ProtoMol::Vector3DBlock velocities,
             ScalarStructure energies,
             Real timestep) {
   self->topology = topo;
   self->positions.vec = positions.vec;
   self->positions.c = positions.c;
   self->velocities.vec = velocities.vec;
   self->velocities.c = velocities.c;
   self->energies = energies;
   self->integrator = new STSIntegrator(timestep, NULL);
   self->outputCache.initialize(self);
}
};
