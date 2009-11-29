#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/config.h>

#include <fah/core/Core.h>
#include <fah/checksum/ChecksumManager.h>
#include <fah/Exception.h>
#include <fah/util/Logger.h>
#include <fah/os/SystemUtilities.h>

#include <iostream>
#include <climits>

using namespace std;
using namespace ProtoMol;
using namespace FAH;

extern void moduleInitFunction(ModuleManager *);

class ProtoMolCore : public Core {
public:
  ProtoMolCore() : Core("ProtoMol", 180, 14) {}

  int init(int argc, char *argv[]) {
    getOptions().add("steps-per-gen")->setType(Option::INTEGER_TYPE);

    int ret = Core::init(argc, argv);
    if (ret) return ret;

    // Validate checkpoint file
    // TODO move this to core.validateCheckpointFile("checkpt");
    if (SystemUtilities::exists("checkpt") &&
        !ChecksumManager::instance().has("checkpt")) {
      // Checksum not valid
      LOG_ERROR("Guru Meditation: checkpt sum");
      return BAD_FRAME_CHECKSUM;
    }

    // Add config file to args
    getArgs().push_front((char *)"protomol.conf");

    return 0;
  }


  int main(int argc, char *argv[]) {
    try {
#if defined(DEBUG) && !defined(_WIN32)
      ProtoMol::Debugger::initStackTrace(argv[0]);
#endif

      ModuleManager modManager;
      moduleInitFunction(&modManager);

      ProtoMolApp app(&modManager);

      app.splash(cout);

      // Load configuration options
      if (!app.load(argc, argv)) return 1;

      // Modify configuration
      // Add outputs FAHFile and FAHGUI
      app.config["FAHFile"] = "../current.xyz";
      app.config["FAHGUI"] = "ProtoMol";

      // Setup checkpointing
      app.config["Checkpoint"] = "checkpt";
      app.config["CheckpointFreq"] = INT_MAX; // Disable

      // Configure and build
      app.configure();
      app.build();

      // Find CheckpointOutput
      OutputCheckpoint *oCheckpt = 0;
      OutputCollection::const_iterator it;
      const OutputCollection *outputs = app.outputs;
      for (it = outputs->begin(); it != outputs->end(); it++)
        if ((oCheckpt = dynamic_cast<OutputCheckpoint *>(*it))) break;

      if (!oCheckpt) THROW("Could not find OutputCheckpoint");

      // Set F@H info
      unsigned gen = getUnit().gen();
      int stepsPerGen = getOptions()["steps-per-gen"].toInteger();
      int firstStep = gen * stepsPerGen;
      int frameSize = stepsPerGen < 100 ? 1 : stepsPerGen / 100;
      setInfo(stepsPerGen, frameSize);

      if (firstStep != app.currentStep) return BAD_WORK_UNIT;

      // Print configuration
      app.print(cout);

      try {
        do {
          // Update shared info file etc.
          step(app.currentStep - firstStep);
          
          if (shouldCheckpoint()) {
            oCheckpt->doIt(app.currentStep);
            checkpoint();
          }
        } while (!shouldQuit() && app.step(min(frameSize, 100)));

      } catch (const ProtoMol::Exception &e) {
        getUnit().type() = CORE_WORK_FAULTY;
      }

      oCheckpt->doIt(app.currentStep);
      app.finalize();
      checkpoint();

    } catch (const ProtoMol::Exception &e) {
      cerr << "ProtoMol ERROR: " << e << endl;
      throw e;
    }

    return 0;
  }
};


int main(int argc, char *argv[]) {
  return ProtoMolCore().run(argc, argv);
}
