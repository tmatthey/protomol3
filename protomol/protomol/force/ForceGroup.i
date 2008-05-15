%module ForceGroup
%{
#include "ForceGroup.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/force/system/SystemForce.h>
using namespace ProtoMol;
%}

%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include "ForceGroup.h"

%extend ProtoMol::ForceGroup {
    void clear() {
      self->getForces().clear();
    }
}