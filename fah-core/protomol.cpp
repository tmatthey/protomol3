#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>

#include <protomol/config.h>
#include <fah/core/main.h>
#include <fah/core/Core.h>
#include <fah/String.h>

#include <iostream>

using namespace std;
using namespace ProtoMol;
using namespace FAH;

extern void moduleInitFunction(ModuleManager *);

extern "C" int core_main(Core &core, const vector<string> &args) {
  try {
#ifdef DEBUG
    ProtoMol::Debugger::initStackTrace(args[0]);
#endif

    ModuleManager modManager;
    moduleInitFunction(&modManager);

    ProtoMolApp app(&modManager);

    app.splash(cout);

    vector<string> moreArgs(args);
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
  return fah_main(argc, argv, core_main);
}
