#include <src/TMDApplication.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/force/ForceGroup.h>

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#endif

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace std;

using std::string;

extern void stringModuleInitFunction(ModuleManager *);

TMDApplication::TMDApplication(const string &tmd_config) {
   tmd_config_filename = tmd_config;
   myApp = NULL;
   fptr = NULL;
}

TMDApplication::~TMDApplication() {

   if (myApp != NULL) {
      myApp->finalize();
      //delete myApp;
   }
}


bool TMDApplication::GetForcePointer(string force_name) {

    if (fptr!=NULL) return false;
    if (myApp == NULL) return false;
    const Integrator *myInt = myApp->integrator;
    while (myInt != NULL) {
       if (IdentifyForce(force_name, myInt) ) {
          report << plain <<"Force pointer obtained from "<<myInt->getId()<<endr;
          return true;
       } else myInt = myInt->next(); // if last integrator, next() returns NULL
    }
    return false;

}

bool TMDApplication::IdentifyForce(string force_name, const Integrator *myInt) {

    ForceGroup *f = myApp->integrator->getForceGroup();
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

void TMDApplication::initialize(string progname, string _fname) {

   //is this okay to do it? Probably yes...
   if (myApp != NULL) return;

   cout<<"Okay until here a..."<<endl;
   cout<<"TMDApplication::tmd_config "<<_fname<<endl;

   //initialize TMD application
   myApp = initialize_TMDApp(progname,_fname);

   //keep force pointer
   string forceName("HarmNMRestForce");

   if (!GetForcePointer(forceName)) {
      cout << "Cant find force pointer in TMDApplication!!"<<endl;
      exit(-1);
   }

   exit(-1);

}

StringProtoMolApp *TMDApplication::initialize_TMDApp(string progname, string _fname) {

   vector<string> myArgs;

   myArgs.push_back(progname);

   if (myApp != NULL) {
       cout << "TMDAPP probably already initialized..."<<endl;
       return NULL;
   }

   ModuleManager m;
   stringModuleInitFunction(&m);
   cout <<"Good until here a "<<endl;
   StringProtoMolApp *app = new StringProtoMolApp(&m);
   cout <<"Good until here b "<<endl;
   cout<<myArgs[0]<<endl;
   myArgs.push_back(_fname);
   cout <<"Good until here ca "<<endl;
   app->configure(myArgs);

   cout <<"Good until here c "<<endl;
   app->build();
   cout <<"Good until here d "<<endl;

   return app;

}


//
//This function will assign a position vector (Vector3DBlock) to myApp
//
void TMDApplication::AssignPos(const Vector3DBlock &pos) {

  myApp->positions.intoAssign(pos);

}

//
//This function will re-initialize spring constants to kappaScale and
//assign target C values
//
void TMDApplication::Reset_target_C(double *cc) {

   fptr->SetCHat(cc,1);

   fptr->SetEigenvalues();

}


//
//Run simulations
//
void TMDApplication::run() {

   while (myApp->step()) continue;

}
