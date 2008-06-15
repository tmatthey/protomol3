%module PySystemForce
%{
#include <protomol/force/system/PySystemForce.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include <protomol/force/system/PySystemForce.h>
