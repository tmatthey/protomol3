#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/base/Exception.h>
#include <protomol/output/OutputCollection.h>
#include <protomol/io/DCDTrajectoryReader.h>
//#include <protomol/config/InputValue.h>
#include <protomol/base/Report.h>

#include <iostream>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

extern void moduleInitFunction(ModuleManager *);

//namespace ProtoMol{
//  declareInputValue(InputDcd, STRING, NOTEMPTY);
//  defineInputValue(InputDcd, "inputDcdFile");
//}

int main(int argc, char *argv[]) {  

  if(argc < 3){
    report << plain << "Only " << argc - 1 << " input parameter, require .conf and input dcd file names." << endr;
    return 0;
  }
  string file_name( argv[2] );
  // dcd file
  DCDTrajectoryReader in;
  try {
    TimerStatistic::timer[TimerStatistic::WALL].start();
    ModuleManager modManager;
    moduleInitFunction(&modManager);
    ProtoMolApp app(&modManager);

    vector<string> args;
    args.push_back("ProtoMol");
    args.push_back(argv[1]);
    if (!app.configure(args)) return 0;

    app.splash(cout);
    app.build();
    if ((int)app.config[InputDebug::keyword]) app.print(cout);

    //open dcd file
    if(in.open(file_name)){ 
      try{
        in >> app.positions;
      }catch(const Exception &e){
        app.currentStep = app.lastStep;
      }
    }else{    
      report << plain << "Input file not defined." << endr;
      return 0;
    }
    
    while(app.currentStep < app.lastStep){
      app.outputs->run(app.currentStep);

      int inc = app.outputs->getNext() - app.currentStep;
      inc = std::min(app.lastStep, app.currentStep + inc) - app.currentStep;
      app.currentStep += inc;

      try{
        in >> app.positions;
      }catch(const Exception &e){
        app.currentStep = app.lastStep;
      }

    }
    //
    app.finalize();
    TimerStatistic::timer[TimerStatistic::WALL].stop();
    return 0;

  } catch (const Exception &e) {
    cerr << "ERROR: " << e << endl;
  }

  return 1;
}



