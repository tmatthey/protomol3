
#include <protomol/base/ModuleManager.h>
#include <protomol/output/OutputCollection.h>
#include <src/integrator/normal/NormalModeStringDiag.h>
#include <protomol/topology/TopologyUtilities.h>

#include <src/StringMethod.h>

//#define DEBUG_METHOD


#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#endif

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace std;

extern void stringModuleInitFunction(ModuleManager *);

StringMethod::StringMethod() : myReparameterizer() {
   myPhiPsiOutput = NULL;
   ppp = NULL;
   eigvectors.clear();

   tmd = NULL;

}

/* // Think about how should I construct these objects */
StringMethod::~StringMethod() {
   if (ppp != NULL) delete [] ppp;
   eigvectors.clear();

   if (tmd != NULL) tmd->finalize();

}

//
// Initialization routine for string simulation. 
// filename : input file containing names of the config files for each 
//            intermediate point of the string.
//
// thres_hold : to find a common space.
//
void StringMethod::initialize(string filename, string progname, double thres_hold, int do_endpoints, string tmd_config) {

    _progname = progname;
    _tmd_confname = tmd_config;

    if (!configReader.open(filename)) report<<error<<"Cant open string config file"<<endr;
 
    if (!configReader.read()) report<<error<<"Error reading string config file"<<endr;

    numpoints = configReader.configfiles.size();

   //Get pointer to config filenames (strings)
   vector<string> *filenames = configReader.GetConfigFilenames();
    
   initializeString(progname, filenames);
   _N = apps[0]->positions.size();
   _3N = apps[0]->positions.size()*3;


   NormalModeStringDiag *nmInt = myIntegratorPointers[0];

   //initialize NormalModeProjection
   //for projection in the subspace defined by eigenvectors/
   //subspace coordinates. 
   
   normalModeProjection.numEigvectsu = _3N;
   cout <<" StringMethod initialize : first and num mode to normalMode Projection "<<nmInt->firstMode<<" "<<nmInt->numMode<<endl;
   normalModeProjection.initialize_utilities(nmInt->firstMode,nmInt->numMode, nmInt->myGamma, nmInt->mySeed, nmInt->myTemp);
   Vector3DBlock *myF = nmInt->getForces();
   normalModeProjection.initialize(_N,apps[0]->topology, myF, true);

   dims = nmInt->GetDimension();

   double *eigvec0 = nmInt->mhQu;

   report << plain <<"StringMethod::initialize : dims "<<dims<<", numpoints "<<numpoints<<", _3N "<<_3N<<endl;

   pos0.intoAssign(apps[0]->positions);
   
   //
   //Initialize CommonSpace object
   //

   if (thres_hold != 0) {
      myCommonSpace.initialize(numpoints, _3N, dims, thres_hold, eigvec0);

      //initialize eigenvector pointers
      initialize_eigenvector_pointers();

      //Compute common basis set
      myCommonSpace.run(eigvectors);

      basis_dim = myCommonSpace.basis_size;
   } else {
      basis_dim = dims;
   }
   int reparam_dim = nmInt->numMode;

   myReparameterizer.initialize(tmd_config, numpoints, reparam_dim, _3N);

   tmd = initialize_StringProtoMolApp(progname, tmd_config, 1);

   double *ev;
   if (thres_hold != 0) 
      ev = myCommonSpace.commonBasis;
   else  //use eigenvectors
      ev = eigvec0;

   normalModeProjection.CopyCommonBasis(ev, basis_dim);

   //calculate C values
   Get_C_Coordinates();

   myReparameterizer.runReparam();

   FixPositionAfterReparameterize();

   if (do_endpoints) {
      start_point = 0;
      end_point = numpoints;
   }else {
       start_point = 1;
       end_point = numpoints-1;
   }


   myPhiPsiOutput = new StringPhiPsiOutput("sp.out",numpoints);
   myPhiPsiOutput->open();

   if (ppp == NULL) ppp = new double[numpoints*2];

   GetPhiPsiCoordinates(); //finds Phi-Psi angles.
   myPhiPsiOutput->run(ppp);

}

void StringMethod::GetCommonBasis() {

   eigvectors.clear();
   initialize_eigenvector_pointers();

   myCommonSpace.run(eigvectors);

   basis_dim = myCommonSpace.basis_size;
}

bool StringMethod::StringStep() {

#if 0
   for(int kk=0;kk<numpoints;kk++) {
      cout<<"apps["<<kk<<"]->positions.size() = "<<apps[kk]->positions.size()<<endl;
   }
#endif

   int currentStep = apps[0]->currentStep;

   int lastStep = apps[0]->lastStep;

   if (currentStep >= lastStep) return false;

   OutputCollection *o;

   for(int i=0;i<numpoints;i++) {
      o = apps[i]->outputs;
      o->run(currentStep);
   }

#if 0
   for(int kk=0;kk<numpoints;kk++) {
      cout<<"apps["<<kk<<"]->positions.size() aa = "<<apps[kk]->positions.size()<<endl;
   }
#endif
   
   o = apps[0]->outputs;

   //calculate inc increment of steps from last o pointer
   int inc = o->getNext() - currentStep;

   //cout<<"StringStep : inc "<<inc<<endl;

   inc = std::min(lastStep, currentStep + inc) - currentStep;


   int p = (int)apps.size();

   for(int i=0;i<p;i++)   apps[i]->currentStep += inc;

   NormalModeStringDiag *nm;

   for(int k=0;k<inc;k++) {

      for(int i=start_point;i<end_point;i++) {
         cout<<"************************* Simulating point "<<i<<"******************"<<endl;
         nm = myIntegratorPointers[i];
         nm->run(1);

      }

      //code for common space might go in here

      Get_C_Coordinates();

      cout<<"*********** Reparameterizing ****************************************"<<endl;
      myReparameterizer.runReparam();

      FixPositionAfterReparameterize();
      
   }

   return true;

}

//
// filenames : config filenames for intermediate points
//
//  This function creates a StringProtoMolApp object for each intermediate point of the string.
//  The array of StringProtoMolApp objects are kept in apps. This is a serial implementation of
//  string method. The idea is the following : during each simulation step, the integrator
//  of each apps object will be executed. I keep a pointer to the set of integrators in
//  myIntegratorPointers.
//

void StringMethod::initializeString(string progname, vector<string> *filenames) {

   unsigned int sz = filenames->size();

   for (unsigned int i=0;i<sz;i++) {
       string s((*filenames)[i]);
       StringProtoMolApp * a = initialize_StringProtoMolApp(progname, s, 0);
       apps.push_back(a);
   }
}

StringProtoMolApp *StringMethod::initialize_StringProtoMolApp(string pname, string filename, int tmd_app) {

   vector<string> myArgs;

   myArgs.push_back(pname);
   myArgs.push_back(filename);
 
   ModuleManager m;
   stringModuleInitFunction(&m);

   StringProtoMolApp *app = new StringProtoMolApp(&m);

   app->configure(myArgs);

   app->build();
   
   if (!tmd_app) {
      NormalModeStringDiag *nmInt = dynamic_cast<NormalModeStringDiag *>(app->integrator);
      if (nmInt == NULL) report <<error<<"No normalmodestringupdate integrator in ProtoMolApp "<<endr;
      myIntegratorPointers.push_back(nmInt);
   }else {
      app->SetFirstStep();
      string myForce("HarmNMRestForce");
      app->GetForcePointer(myForce);
   }
   myArgs.pop_back();
   myArgs.pop_back();

   return app;

}


//
//Only needed for common vectors.
//

void StringMethod::initialize_eigenvector_pointers() {

   double *ev1 = myIntegratorPointers[0]->mhQu;
   unsigned int sz = myIntegratorPointers.size();
   double *ev2 = myIntegratorPointers[sz-1]->mhQu;
   eigvectors.push_back(ev1);
   eigvectors.push_back(ev2);
}

void StringMethod::OutEigVec() {

   for (unsigned int i=0;i<myIntegratorPointers.size();i++) {
      //double *evec = myIntegratorPointers[i]->mhQu;

   }
}

void StringMethod::Get_C_Coordinates() {

   v0.intoAssign(pos0);

   for(int i=0;i<numpoints;i++) {
      vi.intoAssign(apps[i]->positions);

      double *cc = myReparameterizer.SubspaceCoordinate(i);

      normalModeProjection.doModeSpaceProj(cc, &vi, &v0, i);
      
   }

}


void StringMethod::FixPositionAfterReparameterize() {


    cout<<"starting FixPositionAfterReparameterize()"<<endl;
    //cout<<"v0.size() "<<v0.size()<<", apps[0]->positions.size() "<<apps[0]->positions.size()<<", apps[0] "<<apps[0]<<" "<<&(apps[0]->positions)<<endl;

    //v0.intoAssign(apps[0]->positions);
    v0.intoAssign(pos0);
    //cout<<"starting FixPositionAfterReparameterize() a"<<endl;

    //cout<<"Zerov.size() "<<ZeroV.size()<<", v0.size() "<<v0.size()<<endl;
    ZeroV.resize(v0.size());
    ZeroV.zero();

    for(int i=1;i<numpoints-1;i++) {
        vi.intoAssign(apps[i]->positions);
        //vi.intoAssign(apps[0]->positions);

        xbar.intoAssign(apps[i]->positions);
        xbar.intoSubtract(v0); //x-x_0
        normalModeProjection.nonSubspacePosition(&xbar, &xbar); //got \bar{x}
       
        //fixing \hat{x}
        double *c_hat = myReparameterizer.GetS(i);
        normalModeProjection.cartSpaceProj(c_hat, &vi, &ZeroV); //vi has modified \hat{x} after reparameterization
        vi.intoAdd(xbar);
        vi.intoAdd(v0);

        //
        //Now...call the TMD method
        //
        //myReparameterizer.SetTMDApp(_progname, vi, i, _tmd_confname);
        
        //
        // Run targeted MD
        //
        //myReparameterizer.RunTMD();

        runTMD(vi, c_hat);

        apps[i]->positions.intoAssign(vi);
         
    }

}

void StringMethod::GetPhiPsiCoordinates() {

  for(int i=0;i<numpoints;i++) {
     ppp[2*i+0] = computePhiDihedral(apps[i]->topology,&(apps[i]->positions),10);
     ppp[2*i+1] = computePhiDihedral(apps[i]->topology,&(apps[i]->positions),17);
  }

}

void StringMethod::ResetTMD() {

   tmd->ResetTime();

}

void StringMethod::UpdateTMDPositionsAndTarget(const Vector3DBlock &v_vec, double *cc) {

   if (tmd->fptr == NULL) {
      cout<<"Error : pointer to harmonic restraint not available."<<endl;
      exit(-1);
   }

   tmd->positions.intoAssign(v_vec); /* v is the position of an intermediate point */

   tmd->fptr->SetCHat(cc,1);

   tmd->fptr->SetEigenvalues(); 
}

void StringMethod::runTMD(const Vector3DBlock &v_vec, double *cc) {

   ResetTMD();

   UpdateTMDPositionsAndTarget(v_vec, cc);

   while (tmd->step()) continue;

}    
