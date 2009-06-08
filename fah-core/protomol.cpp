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

#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
using namespace ProtoMol;
using namespace FAH;

extern void moduleInitFunction(ModuleManager *);

extern "C" int core_main(int argc, char *argv[]) {
  try {
#ifdef DEBUG
    ProtoMol::Debugger::initStackTrace(argv[0]);
#endif

    FAH::Core &core = *FAH::Core::getInstance();

    ModuleManager modManager;
    moduleInitFunction(&modManager);

    ProtoMolApp app(&modManager);

    app.splash(cout);

    // Add outputs FAHFile and FAHGUI
    string name = "ProtoMol Project";
    name += core.unit->project_id;
    app.config["FAHFile"] = "../current.xyz";
    app.config["FAHGUI"] = name;

    // Setup checkpointing
    app.config["Checkpoint"] = "checkpt";
    app.config["CheckpointFreq"] = -1; // Disable


    if (!app.configure(argc, argv)) return 1;
    app.build();

    // Find CheckpointOutput
    OutputCheckpoint *oCheckpt = 0;
    OutputCollection::const_iterator it;
    const OutputCollection *outputs = app.outputs;
    for (it = outputs->begin(); it != outputs->end(); it++)
      if ((oCheckpt = dynamic_cast<OutputCheckpoint *>(*it))) break;

    if (!oCheckpt) THROW("Could not find OutputCheckpoint");

    app.print(cout);


    // FAH Core setup
    // Initialize shared file
    core.initSharedInfo(name, app.lastStep, 1);

    int outputFreq = toInt(app.config["outputfreq"]);
    while (app.step()) {
      if (app.currentStep % outputFreq == 0)
        cout << "Step: " << app.currentStep << endl;

      // Update shared info file.
      core.updateSharedInfo(app.currentStep);

      if (core.doCheckpoint()) {
        oCheckpt->doIt(app.currentStep);
        core.checkpoint();
      }
    }
      
    oCheckpt->doIt(app.currentStep);
    core.checkpoint();
    app.finalize();

    // Return results
    core.addResultFiles("*");
    core.addResultFile(String::printf("../logfile_%s.txt",
                                      core.suffix.c_str()));

    // Setup header
    core.unit->core_type = 180;
    core.unit->compression_type = 1;

    return 0;

  } catch (const ProtoMol::Exception &e) {
    cerr << "ProtoMol ERROR: " << e.getMessage() << endl;
#ifdef DEBUG
    cerr << e << endl;
#endif
  }

  return 1;
}

int main(int argc, char *argv[]) {
  int ret;

  try {
    Core core;

    ret = core.init(argc, argv);
    if (ret) return ret;

    // Validate checkpoint file
    struct stat buf;
    if (!stat("checkpt", &buf) && !ChecksumManager::instance().has("checkpt")) {
      // Checksum not valid
      cerr << "Guru Meditation: checkpt sum" << endl;
      return BAD_FRAME_CHECKSUM;
    }

    // Add config file
    core.args.push_back((char *)"protomol.conf");

    if (core_main(core.args.size(), &core.args[0]))
      return UNKNOWN_ERROR;

    ret = core.finalize();

  } catch (const FAH::Exception &e) {
    cerr << "Core ERROR: " << e << endl;

    if (e.getCode() == BAD_ARGUMENTS) Core::usage(argv[0]);

    if (e.getCode()) return e.getCode();
    else return UNKNOWN_ERROR;

  } catch (const ProtoMol::Exception &e) {
    cerr << "Core ERROR: " << e << endl;

    return UNKNOWN_ERROR;
  }

  return ret;
}
