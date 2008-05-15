%module AngleForce
%{
#include <protomol/topology/Angle.h>
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include "AngleSystemForce.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include "AngleSystemForce.h"
%template(ASF_PBC) ProtoMol::AngleSystemForce<ProtoMol::PeriodicBoundaryConditions>;
%template(ASF_VBC) ProtoMol::AngleSystemForce<ProtoMol::VacuumBoundaryConditions>;