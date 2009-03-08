#include <src/Reparameterizer.h>
#include <iostream>


using namespace ProtoMol;
using namespace std;

Reparameterizer::Reparameterizer(int np, int r_d, int _3n) :
  numpoints(np), reparam_dim(r_d), _3N(_3n) {
   //tmdApp = NULL;

   myReparam = NULL;
}

Reparameterizer::~Reparameterizer() {

   //if (tmdApp != NULL) delete tmdApp;

   if (myReparam != NULL) delete myReparam;
}

//
//tmd_config is the name of the config file for TMD application.
//I only need this when I am initializing TMD application. Thats why it ok
//to just pass this to the initialization routine. I dont need to keep this
//as a member variable.
//

//
//This function initializes the classes required for reparameterization
//
void Reparameterizer::initialize(string tmd_config, int np, int r_d, int _3n) {

   numpoints = np;
   reparam_dim = r_d;
   _3N = _3n;

   myReparam = new Reparam(numpoints, reparam_dim);

   //initialize reparam object...just once here
   myReparam->initialize(_3N);

   myReparam->initialize_stringPos_DataStructures(numpoints, _3N);

   cout<<"Reparameterizer:: tmd_config "<<tmd_config<<endl;

   //if (tmdApp == NULL) tmdApp = new TMDApplication(tmd_config);

}

//
//Returns the array pointer which holds subspace coordinate from simulation for
//intermediate point i. Once these arrays are filled on, run reparameterization
//(piecewiese linear interpolation + targeted MD)
//
double *Reparameterizer::SubspaceCoordinate(int i) {

  return myReparam->stringPos[i];

}

//
// Set up TMD application to move to the next target point in the subspace.
//
void Reparameterizer::SetTMDApp(string _progname, const Vector3DBlock &v, int i, string _fname) {
   
  //tmdApp->initialize(_progname, _fname);

  //tmdApp->AssignPos(v);

  //tmdApp->Reset_target_C(myReparam->S[i]);

}

//
//Run piecewise linear interpolation.
//
void Reparameterizer::runReparam() {

   myReparam->run();

}

//
//Run targeted MD.
//
void Reparameterizer::RunTMD() {

  //tmdApp->run();

}
