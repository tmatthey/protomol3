%module DCDTrajectoryReader
%{
#include "DCDTrajectoryReader.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include "DCDTrajectoryReader.h"
