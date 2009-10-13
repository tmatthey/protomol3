#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/output/OutputCheckpoint.h>
#include <protomol/config.h>

#include <fah/core/Core.h>
#include <fah/core/chksum/ChecksumManager.h>
#include <fah/Exception.h>
#include <fah/String.h>
#include <fah/util/Logger.h>
#include <fah/os/SystemUtilities.h>

#include <iostream>
#include <climits>

using namespace std;
using namespace ProtoMol;
using namespace FAH;

extern void moduleInitFunction(ModuleManager *);

extern "C" int core_main(int argc, char *argv[]) {
  try {
#if defined(DEBUG) && !defined(_WIN32)
    ProtoMol::Debugger::initStackTrace(argv[0]);
#endif

    FAH::Core &core = *FAH::Core::getInstance();

    ModuleManager modManager;
    moduleInitFunction(&modManager);

    ProtoMolApp app(&modManager);

    app.splash(*LOG_RAW_STREAM());

    // Load configuration options
    if (!app.load(argc, argv)) return 1;

    // Modify configuration
    // Add outputs FAHFile and FAHGUI
    string name = "ProtoMol Project";
    name += core.unit->project_id;
    app.config["FAHFile"] = "../current.xyz";
    app.config["FAHGUI"] = name;

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

    // Print configuration
    app.print(*LOG_RAW_STREAM());

    // FAH Core setup
    // Initialize shared file
    core.initSharedInfo(name, app.lastStep, 1);

    int outputFreq = toInt(app.config["outputfreq"]);
    while (!core.shouldQuit() && app.step(100)) {
      if (app.currentStep % outputFreq == 0)
        LOG_INFO(1, "Step: " << app.currentStep);

      // Update shared info file.
      core.updateSharedInfo(app.currentStep);

      if (core.doCheckpoint()) {
        oCheckpt->doIt(app.currentStep);
        core.checkpoint();
      }
    }
      
    oCheckpt->doIt(app.currentStep);
    app.finalize();
    core.checkpoint();

    // Return results
    core.addResultFiles("*");

    // Setup header
    core.unit->core_type = 180;
    core.unit->compression_type = 1;

    return 0;

  } catch (const ProtoMol::Exception &e) {
    LOG_ERROR("ProtoMol ERROR: " << e);
  }

  return 1;
}

int main(int argc, char *argv[]) {
  int ret;

  try {
    Core &core = *FAH::Core::getInstance();

    ret = core.init(argc, argv);
    if (ret) return ret;

    // Validate checkpoint file
    if (SystemUtilities::exists("checkpt") &&
        !ChecksumManager::instance().has("checkpt")) {
      // Checksum not valid
      LOG_ERROR("Guru Meditation: checkpt sum");
      return BAD_FRAME_CHECKSUM;
    }

    // Add config file
    vector<char *> args;
    vector<char *>::iterator it = core.args.begin();
    args.push_back(*it++); // Executable name
    args.push_back((char *)"protomol.conf");
    args.insert(args.end(), it, core.args.end());
    args.push_back(0); // Sentinel

    if (core_main(args.size() - 1, &args[0])) return UNKNOWN_ERROR;

    ret = core.finalize();

  } catch (const FAH::Exception &e) {
    if (e.getCode() == BAD_ARGUMENTS)
      Core::getInstance()->usage(cerr, argv[0]);

    LOG_ERROR("Core: " << e);

    if (e.getCode()) return e.getCode();
    else return UNKNOWN_ERROR;

  } catch (const ProtoMol::Exception &e) {
    LOG_ERROR("std::exception: " << e);

#ifdef DEBUG
    throw e; // Rethrow to get core dump
#endif

    return UNKNOWN_ERROR;
  }

  return ret;
}
