%module PDBReader
%{
#include "PDBReader.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include "std_string.i"
//%include "std_vector.i"
%include "../type/Vector3DBlock.h"
%include "../type/PDB.h"
%include "File.h"
%include "Reader.h"
%include "PDBReader.h"

