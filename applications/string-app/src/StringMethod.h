/*  -*- c++ -*-  */
#ifndef STRINGMETHOD_H
#define STRINGMETHOD_H


#include <src/StringProtoMolApp.h>
#include <src/integrator/normal/NormalModeStringDiag.h>

#include <src/base/StringConfigReader.h>
#include <src/base/CommonSpace.h>
#include <src/Reparameterizer.h>
#include <src/NormalModeProjection.h>
#include <src/base/StringPhiPsiOutput.h>
#include <vector>

class StringMethod {

   public:

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //  constructors and destructors
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      StringMethod();

      ~StringMethod();


      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Member functions
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public :

      void initialize(string filename, string progname, double thres_hold, int do_endpoints, string tmd_config);

      bool StringStep();
      void GetPhiPsiCoordinates();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Private member functions
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   private:
       
      void initializeString(string progname, vector<string> *filenames);
   
      StringProtoMolApp *initialize_StringProtoMolApp(string progname, string filename, int tmd_app);

      void initialize_eigenvector_pointers();

      void Get_C_Coordinates();

      void GetCommonBasis();

      void FixPositionAfterReparameterize();

      void OutEigVec();

   public:
 
      //
      // Methods for manipulating tmd object
      //

      void ResetTMD();

      void UpdateTMDPositionsAndTarget(const Vector3DBlock &v_vec, double *cc);

      void runTMD(const Vector3DBlock &v_vec, double *cc);

   public:

      std::vector<StringProtoMolApp *> apps; /* ProtoMol applications...each string point is an app */

      StringConfigReader configReader; /* string config file reader (contains a list of configs
                                          for intermediate points */

      int numpoints ; /* no of intermediate points on the string */

      int dims ; /* subspace size at each intermediate point */

      std::vector<NormalModeStringDiag *> myIntegratorPointers; /* pointers to integrator objects at
                                                             each intermediate point */

      CommonSpace myCommonSpace;   

      Reparameterizer myReparameterizer;

      StringPhiPsiOutput *myPhiPsiOutput;
   
      double *ppp;
      std::vector<double *> eigvectors;


      int basis_dim;

      int _N, _3N;
                                      
      NormalModeProjection normalModeProjection;

      Vector3DBlock v0, vi, xbar, ZeroV;

      Vector3DBlock pos0;

      int start_point, end_point;

   private:
      std::string _progname, _tmd_confname;

   public:
      StringProtoMolApp *tmd;
};

#endif /* STRINGMETHOD_H */
