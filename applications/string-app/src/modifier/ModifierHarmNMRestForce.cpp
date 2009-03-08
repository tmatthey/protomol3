#include <src/modifier/ModifierHarmNMRestForce.h>
#include <protomol/topology/Topology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/force/ForceGroup.h>

using namespace ProtoMol::Report;
using std::deque;
using std::vector;
using std::string;

namespace ProtoMol {

    //  _________________________________________________________ ModifierHarmNMRestForce

    ModifierHarmNMRestForce::ModifierHarmNMRestForce(const Integrator *i) : Modifier(),
          myIntegrator(i){
       fptr = NULL;
    }

    void ModifierHarmNMRestForce::doInitialize() {

        report << plain <<"Initializing HarmNM Modifier"<<endr;

       //find the pointer to HarmNMRestForce
        string forceName("HarmNMRestForce");
           
       if (!GetForcePointer(forceName)) report << error <<"Cant find HarmNMForce "<<endr; 

       enable();

    }

    Real ModifierHarmNMRestForce::getTimestep() const {
       return (myIntegrator->getTimestep());
    }

    bool ModifierHarmNMRestForce::GetForcePointer(string force_name) {

        if (fptr!=NULL) return false;
        const Integrator *myInt = myIntegrator;
        while (myIntegrator != NULL) {
           if (IdentifyForce(force_name, myInt) ) {
              report << plain <<"Force pointer obtained from "<<myInt->getId()<<endr;
              return true;
           } else myInt = myInt->next(); // if last integrator, next() returns NULL
        }
       return false; 
    }

    bool ModifierHarmNMRestForce::IdentifyForce(string force_name, const Integrator *myInt) {

        ForceGroup *f = myIntegrator->getForceGroup();
        vector<Force *> myForces = f->getForces();
        HarmNMForce *p;
        for(unsigned int i=0 ;i<myForces.size();i++) {
            p = dynamic_cast<HarmNMForce *>(myForces[i]);
            if (p) {
                fptr = p;
                return true;
            }
        }
        return false;
    }

    void ModifierHarmNMRestForce::doExecute(Integrator *i) {

       //report << " Running do Execute of HarmNM modifier"<<endr;

       int step = (int)(app->topology->time/getTimestep());
       fptr->SetStep(step);

    }

}
