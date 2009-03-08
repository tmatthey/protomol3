
#include <src/base/HarmonicForceUtilities.h>
#include <protomol/force/ForceGroup.h>

using namespace std;
using namespace ProtoMol::Report;

namespace ProtoMol {

HarmNMForce *GetHarmNMForcePointer(std::string force_name, const Integrator *integrator) {

   HarmNMForce *fptr = NULL;

   const Integrator *myInt = integrator;
   while ((myInt != NULL) && (fptr == NULL)) {

      fptr = IdentifyHarmNMForce(force_name, myInt);

      myInt = myInt->next();
   }
      
   return fptr;
}

HarmNMForce *IdentifyHarmNMForce(std::string force_name, const Integrator *integrator) {

   ForceGroup *f = integrator->getForceGroup();

   vector<Force *> myForces = f->getForces();

   HarmNMForce *p = NULL;

   for(unsigned int i=0 ;i<myForces.size();i++) {
      p = dynamic_cast<HarmNMForce *>(myForces[i]);
      if (p != NULL) return p;
    }

    return p;

}

}
