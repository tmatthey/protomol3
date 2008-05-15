%module PARReader
%{
#include "PARReader.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include "../type/PAR.h"
%include "PARReader.h"
