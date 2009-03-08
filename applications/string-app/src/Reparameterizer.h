#ifndef REPARAMETERIZER_H
#define REPARAMETERIZER_H

#include <src/Reparam.h>
#include <protomol/type/Vector3DBlock.h>

#include <string>


class Reparameterizer {

  public: 

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //   Constructors, desctructors
     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Reparameterizer() {}
     Reparameterizer(int np, int r_d, int _3n);

     ~Reparameterizer();

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     // Member functions 
     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     void initialize(std::string tmd_conf, int np, int r_d, int _3n);

     double *SubspaceCoordinate(int i);


     double *GetS(int i) { return myReparam->S[i]; }

     void SetTMDApp(std::string _progname, const ProtoMol::Vector3DBlock &v, int i, std::string _fname);

     void runReparam();

     void RunTMD();

  private:

    int numpoints, reparam_dim, _3N;

    //ProtoMol::TMDApplication *tmdApp;
  
    Reparam *myReparam;


};


#endif
