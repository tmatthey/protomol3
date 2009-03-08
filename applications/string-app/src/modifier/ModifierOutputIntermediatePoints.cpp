#include <src/modifier/ModifierOutputIntermediatePoints.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/Topology.h>


using namespace ProtoMol::Report;
using std::deque;
using std::vector;
using std::string;
 
namespace ProtoMol {

   //  _________________________________________________________ ModifierOutputIntermediateaPoints

   ModifierOutputIntermediatePoints::ModifierOutputIntermediatePoints(const Integrator *i, int n_p) : Modifier(),
     myIntegrator(i), next_point(n_p) {
       fptr = NULL;
       //pdbWriter = NULL;
    }

    ModifierOutputIntermediatePoints::~ModifierOutputIntermediatePoints() {
       //if (pdbWriter != NULL) delete pdbWriter;
    }

    void ModifierOutputIntermediatePoints::doInitialize() {

       report << plain <<"Initializing OutputIntermediatePoints Modifier"<<endr;

       //find the pointer to HarmNMRestForce
        string forceName("HarmNMRestForce");


       if (!GetForcePointer(forceName)) report << error <<"Cant find HarmNMForce "<<endr;

       //next_point = 1; /* will be initialized through the constructor */

       pts_filename = NextIntermediatePointFilename();
       cout<<"doInitialize::intermediate point filename "<<pts_filename<<endl;
     
    }

    Real ModifierOutputIntermediatePoints::getTimestep() const {
       return (myIntegrator->getTimestep());
    }


    string ModifierOutputIntermediatePoints::NextIntermediatePointFilename() {
       string fname("intermediate_point_");
       fname += toString<int>(next_point);
       fname += ".pdb";
       return fname;
    }

    void ModifierOutputIntermediatePoints::doExecute(Integrator *i) {
     
       //if (fptr->satisfy_maxnorm) the current point has been reached.
       int current_step = (int)(app->topology->time/getTimestep());
       if ((fptr->satisfy_maxnorm) && (current_step)) {
           report << plain <<"Satisfy_maxnorm = 1, running doExecute"<<endr;
           //write the positions into the PDB file.
           //pdbWriter = new PDBWriter(pts_filename);
           //if (!pdbWriter->open()) report << error <<"Cant open intermediate point coordinates file "<<pts_filename<<endr;
           if (!pdbWriter.open(pts_filename)) report << error <<"Cant open intermediate point coordinates file "<<pts_filename<<endr;

           //pdbWriter->setComment("Pdb file for intermediate point "+toString(next_point));
           pdbWriter.setComment("Pdb file for intermediate point "+toString(next_point));

           //copy the minimal image
           const Vector3DBlock *pos = app->outputCache.minimalPositions();


           if (!pdbWriter.write(*pos, app->outputCache.pdb()))
              report << error <<"Cant write intermediate point "<<next_point<<" to coordinate file "<<pts_filename<<endr;

   
           //if all points has been found then terminate the targeted MD.
           if (next_point == fptr->target_point) //should be called target_points
           {
              report <<"Found all the intermediate points"<<endr;
              exit(0);
           }else {
              //else increment next_point and create the new filename 
              //for intermediate point coordinate.
              next_point++;
              pts_filename = NextIntermediatePointFilename();
              

              //set next target
              fptr->Set_Next_Target();   
              fptr->SetEigenvalues();  
 
              //reset target_steps for the next target
              int this_step = (int)(app->topology->time/getTimestep());
              fptr->next_target_steps = this_step + fptr->target_steps;
           }
        }else {
            //report << plain <<"Max norm condition not satisfied"<<endr;
        }

    }  

    bool ModifierOutputIntermediatePoints::GetForcePointer(string force_name) {
 
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
 
    bool ModifierOutputIntermediatePoints::IdentifyForce(string force_name, const Integrator *myInt) {
 
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


}
