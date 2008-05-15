%module GenericTopology
%{
//#include "CubicCellManager.h"
#include "PeriodicBoundaryConditions.h"
#include "VacuumBoundaryConditions.h"
#include "Topology.h"
#include "BuildTopology.h"
#include <protomol/type/PSF.h>
#include <protomol/type/PAR.h>
using namespace ProtoMol;
%}

%include "../type/Real.h"
%include "std_vector.i"
%include "GenericTopology.h"
%include "SemiGenericTopology.h"
%include "BuildTopology.h"
%template(SGT_Periodic) ProtoMol::SemiGenericTopology <ProtoMol::PeriodicBoundaryConditions >;
%template(SGT_Vacuum) ProtoMol::SemiGenericTopology <ProtoMol::VacuumBoundaryConditions >;


%include "Topology.h"
%template(T_Periodic) ProtoMol::Topology <ProtoMol::PeriodicBoundaryConditions, ProtoMol::CubicCellManager >;
%template(T_Vacuum) ProtoMol::Topology <ProtoMol::VacuumBoundaryConditions, ProtoMol::CubicCellManager >;


%extend ProtoMol::GenericTopology {
        void setExclusion(char* e) {self->exclude = ExclusionType((const char*)e);}

};
%extend ProtoMol::Topology<ProtoMol::VacuumBoundaryConditions, ProtoMol::CubicCellManager> {
	void setCellSize(ProtoMol::Real r){
	   self->cellManager.setCellSize(r);
	}
}

%extend ProtoMol::Topology<ProtoMol::PeriodicBoundaryConditions, ProtoMol::CubicCellManager> {
	void setCellSize(ProtoMol::Real r){
	   self->cellManager.setCellSize(r);
	}

        void setBC(float cB1x, float cB1y, float cB1z,
                   float cB2x, float cB2y, float cB2z,
                   float cB3x, float cB3y, float cB3z,
                   float cOx, float cOy, float cOz) {
           self->boundaryConditions.set(ProtoMol::Vector3D(cB1x, cB1y, cB1z),
                                        ProtoMol::Vector3D(cB2x, cB2y, cB2z),
                                        ProtoMol::Vector3D(cB3x, cB3y, cB3z),
                                        ProtoMol::Vector3D(cOx, cOy, cOz));
        }
};
		
		
	
