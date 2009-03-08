#ifndef TMDAPPLICATION_H
#define TMDAPPLICATION_H

//
// TMDApplication class will keep a StringProtoMolApp object for targeted MD.
// It should be initialized with appropriate config file by StringMethodApp.
// Next, it will initialize the StringProtoMolApp object for TMD and keep a pointer
// to the harmonic restraint force. 

//
// In run method, it should accept current positions of a system and a target 
// C vector (subspace coordinates). It should first set the positions and C vector
// appropriately with the required integrator and force respectively.
// Next, it should reset the integration steps to 0 and run m steps of targeted MD.
//

#include <src/StringProtoMolApp.h>
#include <src/base/HarmonicForceUtilities.h>

#include <string>

namespace ProtoMol {

  class TMDApplication {
 

   public:

      TMDApplication(const std::string &tmd_config);
    
      ~TMDApplication();

      void initialize(std::string progname, std::string _fname);

      void run();

      void AssignPos(const ProtoMol::Vector3DBlock &pos);

      void Reset_target_C(double *cc);

   private:
      ProtoMol::StringProtoMolApp *initialize_TMDApp(std::string s, std::string t);

      bool GetForcePointer(std::string force_name);

      bool IdentifyForce(std::string force_name, const ProtoMol::Integrator *myInt);

   
    public:

       std::string tmd_config_filename;

    public:

       ProtoMol::StringProtoMolApp *myApp;
       ProtoMol::HarmNMForce *fptr;


  };
}

#endif /* TMDAPPLICATION_H */
