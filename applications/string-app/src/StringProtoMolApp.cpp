#include <src/StringProtoMolApp.h>


using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;


StringProtoMolApp::StringProtoMolApp(ModuleManager *modManager):
   ProtoMolApp(modManager) {
   fptr = NULL;
}

StringProtoMolApp::~StringProtoMolApp() {}

void StringProtoMolApp::build1() {
  return;
}

void StringProtoMolApp::build2() {
  return;
}

//firstStep should be set to 0 before this is called
void StringProtoMolApp::ResetTime() {

   currentStep = firstStep;

   topology->time = firstStep * integrator->getTimestep();

}

void StringProtoMolApp::GetForcePointer(string force_name) {

   if (fptr != NULL) report << error << "Force pointer is available."<<endr;

   fptr = GetHarmNMForcePointer(force_name, integrator);

}   
