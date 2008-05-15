%module XYZBinReader
%{
#include "XYZBinReader.h"
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include "std_string.i"
%include "../type/Vector3DBlock.h"
%include "../type/XYZ.h"
%include "File.h"
%include "Reader.h"
%include "XYZBinReader.h"
