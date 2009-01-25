#include <protomol/ProtoMolApp.h>

#include <protomol/base/ModuleManager.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/base/Zap.h>
#include <protomol/base/Report.h>

#include <protomol/module/MainModule.h>
#include <protomol/module/IOModule.h>
#include <protomol/module/ConfigurationModule.h>

#include <protomol/type/String.h>

#include <protomol/config/CommandLine.h>
#include <protomol/config/Configuration.h>

#include <protomol/io/ConfigurationReader.h>

#include <protomol/factory/TopologyFactory.h>
#include <protomol/factory/OutputFactory.h>

#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/BuildTopology.h>
#include <protomol/topology/TopologyUtilities.h>

#include <protomol/output/OutputCollection.h>

#include <iomanip>
#ifdef HAVE_PACKAGE_H
#include <protomol/package.h>
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;

ProtoMolApp::ProtoMolApp(ModuleManager *modManager) :
  modManager(modManager), cmdLine(&config), outputs(0), integrator(0),
  topology(0) {
  modManager->init(this);

  topologyFactory.registerAllExemplarsConfiguration(&config);
  outputFactory.registerAllExemplarsConfiguration(&config);
}

ProtoMolApp::~ProtoMolApp() {}

void ProtoMolApp::splash(ostream &stream) {
  const int w = 16;
  stream
    << headerRow("ProtoMol") << endl
    << setw(w) << "Description: ";
  fillFormat(stream, "A rapid PROTOtyping MOLecular dynamics object-oriented "
             "component based framework.", w, w);
  stream 
#ifdef HAVE_PACKAGE_H
    << setw(w) << "Version: " << PACKAGE_VERSION << endl
    << setw(w) << "SVN revision: " << PACKAGE_REVISION << endl
    << setw(w) << "Repository: " << PACKAGE_SOURCE_REPO << endl
    << setw(w) << "Homepage: " << PACKAGE_HOMEPAGE << endl
    << setw(w) << "Report bugs to: " << PACKAGE_BUGREPORT << endl
    << setw(w) << "Compiler: " << PACKAGE_COMPILER << " "
    << PACKAGE_COMPILER_VERSION << endl
    << setw(w) << "Flags: " << PACKAGE_COMPILER_FLAGS << endl
    << setw(w) << "Extra libs: " << PACKAGE_COMPILER_LIBS << endl
    << setw(w) << "Built by: " << PACKAGE_BUILT_BY << endl
    << setw(w) << "Build platform: " << PACKAGE_PLATFORM << endl
#endif // HAVE_PACKAGE_H
    << setw(w) << "Build date: " <<  __DATE__ << ", " << __TIME__ << endl
#ifdef BUILD_FOR_FAH
    << setw(w) << "Checksumming: "
    << "Enabled for Folding@Home file protection." << endl
#endif // BUILD_FOR_FAH
#ifdef HAVE_PACKAGE_H
    << setw(w) << "Please cite: ";
  fillFormat(stream, PACKAGE_CITE, w, w);
  stream
#endif // HAVE_PACKAGE_H
    << PROTOMOL_HR << endl;
}

void ProtoMolApp::configure(const string &configfile) {
  vector<string> args;

  args.push_back("ProtoMol");
  args.push_back(configfile);
  configure(args);
}


bool ProtoMolApp::configure(int argc, char *argv[]) {
  return configure(vector<string>(argv, argv + argc));
}

bool ProtoMolApp::configure(const vector<string> &args) {
  // Parse command line
  if (cmdLine.parse(args)) return false;

  // Read Config file
  if (config.valid(InputConfig::keyword))
    changeDirectory(config[InputConfig::keyword]);
  else THROW("Configuration file not set.");

  modManager->configure(this);

  return true;
}

void ProtoMolApp::build() {
  // Read data
  modManager->read(this);


  // Build topology
  try {
    topology = topologyFactory.make(&config);
  } catch (const Exception &e) {
    // Try to get some defaults with the postions known ...
    const GenericTopology *prototype =
      topologyFactory.find(config[GenericTopology::keyword].getString());

    if (prototype) {
      vector<Parameter> parameters = prototype->getDefaults(positions);

      for (unsigned int i = 0; i < parameters.size(); i++)
        if (!config.valid(parameters[i].keyword) &&
            parameters[i].value.valid()) {
          config.set(parameters[i].keyword, parameters[i].value);
          report << hint << parameters[i].keyword << " undefined, using "
                 << parameters[i].value.getString() << "." << endr;
        }

      topology = topologyFactory.make(&config);
    }

    if (!topology) throw e;
  }

  // Using SCPISM parameter? Flag or filename
  std:string scpismFile;
  scpismFile.clear();
  if (config[InputDoSCPISM::keyword] || config.valid(InputSCPISM::keyword)) {
    if(config[InputDoSCPISM::keyword])
       topology->doSCPISM = config[InputDoSCPISM::keyword];
    else
       topology->doSCPISM = 4;
    if(config.valid(InputSCPISM::keyword)) scpismFile = (std::string)config[InputSCPISM::keyword];
    report << "SCPISM: doSCPISM set to " << topology->doSCPISM << ", using filename " 
      << config[InputSCPISM::keyword] << "." << endr;

  }

  // Build the topology
  buildTopology(topology, psf, par, config[InputDihedralMultPSF::keyword], scpismFile);


  // Register Forces
  modManager->registerForces(this);


  // Build the integrators and forces
  integrator =
    integratorFactory.make(config[InputIntegrator::keyword], &forceFactory);

  // Setup run paramiters (used for GUI so required here)
  currentStep = config[InputFirststep::keyword];
  lastStep = currentStep + (int)config[InputNumsteps::keyword];

  // Create outputs
  // TODO if !Parallel::iAmMaster() turn off outputs
  if (config[InputOutput::keyword])
    outputs = outputFactory.makeCollection(&config);

  else outputs = new OutputCollection; // Empty collection


  // Post build processing
  modManager->postBuild(this);

  report << plain << "Actual start temperature : "
         << temperature(topology, &velocities) << "K" << endr;


  // Add Integrator Modifiers
  modManager->addModifiers(this);

  // Initialize
  energies.molecularVirial(config[InputMolVirialCalc::keyword]);
  energies.virial(config[InputVirialCalc::keyword]);
  report << plain << "Virial tensor : " << energies.virial() << endr;
  report << plain << "Molecular virial tensor : "
         << energies.molecularVirial() << endr;

  topology->time =
    (Real)config[InputFirststep::keyword] * integrator->getTimestep();

  integrator->initialize(this);
  outputs->initialize(this);
  outputCache.initialize(this);

  // Init cache
  //outputs->addToCache(pdbAtoms); // TODO fix this
  outputCache.add(psf);
  outputCache.add(par);

  // Print Factories
  if ((int)config[InputDebug::keyword] >= 5) {
    cout
      << headerRow("Factories")     << endl
      << headerRow("Configuration") << endl << config            << endl
      << headerRow("Topology")      << endl << topologyFactory   << endl
      << headerRow("Integrator")    << endl << integratorFactory << endl
      << headerRow("Force")         << endl << forceFactory      << endl
      << headerRow("Output")        << endl << outputFactory     << endl;
  }

   // Clear all factories
  topologyFactory.unregisterAllExemplars();
  integratorFactory.unregisterAllExemplars();
  forceFactory.unregisterAllExemplars();
  outputFactory.unregisterAllExemplars();
 
  // Setup run paramiters
  //currentStep = config[InputFirststep::keyword];
  //lastStep = currentStep + (int)config[InputNumsteps::keyword];

  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
}

bool ProtoMolApp::step() {
  if (currentStep >= lastStep) return false;

  outputs->run(currentStep);

  int inc = outputs->getNext() - currentStep;
  inc = std::min(lastStep, currentStep + inc) - currentStep;
  currentStep += inc;

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();

  integrator->run(inc);

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();

  return true;
}

void ProtoMolApp::finalize() {
  outputs->finalize(lastStep);

  // Clean up
  zap(topology);
  zap(integrator);
  zap(outputs);

  report
    << allnodesserial << plain << "Timing: " << TimerStatistic() << "." << endr;
}

void ProtoMolApp::print(ostream &stream) {

  // Output
  stream << headerRow("Outputs") << endl;

  for (OutputCollection::const_iterator itr =
         const_cast<const OutputCollection *>(outputs)->begin();
       itr != const_cast<const OutputCollection*>(outputs)->end(); itr++) {

    stream << "Output " << (*itr)->getId();

    vector<Parameter> parameters;
    (*itr)->getParameters(parameters);
    for (unsigned int i = 0; i < parameters.size(); i++)
      stream << " " << parameters[i].value.getString();
    stream << "." << endl;
  }

  if (!((bool)config[InputOutput::keyword]))
    stream << "All output suppressed!" << endl;


  // Integrator
  stream << headerRow("Integrator") << endl;
  vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
  stream  << InputIntegrator::keyword << " {" << endl;

  for (int i = inter.size() - 1; i >= 0; i--)
    stream << Constant::PRINTINDENT << "Level "
           << i << " " << inter[i].print() << endl;

  stream << "}" << endl;


  // Topology
  stream << headerRow("Topology") << endl;
  stream << topology->print(&positions) << endl;

  stream << PROTOMOL_HR << endl;
}
