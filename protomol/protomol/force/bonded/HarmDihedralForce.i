%module HarmDihedralForce
%{
#include <protomol/topology/SemiGenericTopology.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include "HarmDihedralSystemForce.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/base/Report.h>
%include <protomol/force/Force.h>
%include <protomol/force/system/SystemForce.h>
%include "HarmDihedralSystemForce.h"
%template(HDSF_Periodic) ProtoMol::HarmDihedralSystemForce <ProtoMol::PeriodicBoundaryConditions >;
%template(HDSF_Vacuum) ProtoMol::HarmDihedralSystemForce <ProtoMol::VacuumBoundaryConditions >;