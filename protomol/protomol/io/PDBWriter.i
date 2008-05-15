%module PDBWriter
%{
#include "PDBWriter.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include "std_vector.i"
//%rename ("PDBAtom") PDB::Atom;
%template() std::vector<PDB::Atom>;
%template() std::vector<PDB::Ter>;
%include <protomol/type/Vector3DBlock.i>
//%include <protomol/topology/Atom.h>
//%include <protomol/topology/AtomType.h>
%include <protomol/type/PDB.h>
%include <protomol/io/File.h>
%include <protomol/io/Writer.h>
%include "PDBWriter.h"
