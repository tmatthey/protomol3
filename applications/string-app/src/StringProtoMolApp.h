#ifndef STRINGPROTOMOLAPP_H
#define STRINGPROTOMOLAPP_H

#include <protomol/ProtoMolApp.h>
#include <src/base/HarmonicForceUtilities.h>

namespace ProtoMol {

  class StringProtoMolApp : public ProtoMolApp {

   public:

     StringProtoMolApp() { fptr = NULL; }
     StringProtoMolApp(ModuleManager *modManager);
     ~StringProtoMolApp();

     void build1();

     void build2();

   public:
     //Need to call this function after build() and before step()
     void SetFirstStep() { firstStep = currentStep ; }

     void GetForcePointer(std::string s);

   public:
     void ResetTime();

   private:
     int firstStep;

   public:
     HarmNMForce *fptr;

  };

}

#endif /* STRINGPROTOMOLAPP_H */
