#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/config.h>

#include <fah/core/Core.h>
#include <fah/core/chksum/ChecksumManager.h>
#include <fah/Exception.h>
#include <fah/util/Logger.h>
#include <fah/os/SystemUtilities.h>

#include <iostream>
#include <climits>

using namespace std;
using namespace ProtoMol;
using namespace FAH;

extern void moduleInitFunction(ModuleManager *);

extern "C" int core_main(int argc, char *argv[]) {
  Core &core = Core::getInstance();

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
    unsigned gen = core.getUnit().base()->gen;
    int stepsPerGen = core.getOptions()["steps-per-gen"].toInteger();
    int firstStep = gen * stepsPerGen;
    int frameSize = stepsPerGen < 100 ? 1 : stepsPerGen / 100;
    core.setInfo(stepsPerGen, frameSize);

    // Print configuration
    app.print(cout);

    do {
      // Update shared info file etc.
      core.step(app.currentStep - firstStep);

      if (core.doCheckpoint()) {
        oCheckpt->doIt(app.currentStep - firstStep);
        core.checkpoint();
      }
    } while (!core.shouldQuit() && app.step(min(frameSize, 100)));
      
    oCheckpt->doIt(app.currentStep);
    app.finalize();
    core.checkpoint();

    return 0;

  } catch (const ProtoMol::Exception &e) {
    cerr << "ProtoMol ERROR: " << e << endl;
    core.markFaulty();
  }

  return 1;
}


int main(int argc, char *argv[]) {
  try {
    Core core("ProtoMol", 180);

    core.getOptions().add("steps-per-gen")->setType(Option::INTEGER_TYPE);

    int ret = core.init(argc, argv);
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
    vector<char *> args;
    vector<char *>::const_iterator it = core.getArgs().begin();
    args.push_back(*it++); // Executables name
    args.push_back((char *)"protomol.conf");
    args.insert(args.end(), it, core.getArgs().end());
    args.push_back(0); // Sentinel

    core_main(args.size() - 1, &args[0]);

    return core.finalize();

  } catch (const FAH::Exception &e) {
    LOG_ERROR("Core: " << e);

    if (e.getCode()) return e.getCode();
    return UNKNOWN_ERROR;

  } catch (const std::exception &e) {
    LOG_ERROR("std::exception: " << e.what());

#ifdef DEBUG
    throw e; // Rethrow to get core dump
#endif

    return UNKNOWN_ERROR;
  }
}
