#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>
#include <protomol/config.h>

#include <fah/core/Core.h>
#include <fah/Exception.h>
#include <fah/String.h>

#include <iostream>

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

    vector<string> moreArgs(argv + 1, argv + argc);
    moreArgs.push_back("protomol.conf");
    if (!app.configure(moreArgs)) return 1;

    string name = "ProtoMol Project";
    name += core.unit->project_id;

    // Add outputs FAHFile and FAHGUI
    app.config["FAHFile"] = "../current.xyz";
    app.config["FAHGUI"] = name;

    app.build();
    app.print(cout);

    // FAH Core setup
    // Initialize shared file
    core.initSharedInfo(name, app.lastStep, 1);

    while (app.step()) {
      cout << "Step: " << app.currentStep << endl;

      // Update shared info file.
      core.updateSharedInfo(app.currentStep);

      if (core.doCheckpoint()) {
        std::cout << "Checkpointing..." << std::flush;
        // TODO save checkpoint here
        core.checkpoint();
        std::cout << "done" << std::endl;
      }
    }
      
    app.finalize();

    // Return results
    core.addResultFiles("*");

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

    // TODO checkpointing setup here

    if (core_main(core.args.size() - 1, &core.args[0]))
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
