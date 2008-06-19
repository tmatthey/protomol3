//Updated for standalone GUI
#if defined (HAVE_GUI) || defined (HAVE_LIBFAH)

#include <protomol/output/OutputFAHGUI.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#ifdef HAVE_LIBFAH
#include <fah/core/GUIServer.h>
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

#ifdef HAVE_LIBFAH
using namespace FAH;
#endif

//____ OutputFAHGUI
#ifdef HAVE_LIBFAH
const string OutputFAHGUI::keyword("FAHGUI");
#else
const string OutputFAHGUI::keyword("Gui");
#endif

OutputFAHGUI::OutputFAHGUI() : name("ProtoMol"), server() {}

OutputFAHGUI::OutputFAHGUI(const string &name, int freq, int port, int prange, const string &projn) :
  Output(freq), name(name), myPort(port), myPortRange(prange), myProjName(projn), server(0) {}

void OutputFAHGUI::doInitialize() {
#ifdef HAVE_LIBFAH
  server = new GUIServer(myProjName.c_str(), app->topology->atoms.size(),
                          app->topology->bonds.size());
#else
  server = new GUIServer(myProjName.c_str(), app->topology->atoms.size(),
                          app->topology->bonds.size(),myPort,myPortRange);
  server->info.iterations = app->lastStep / 1000;
  server->info.frames = app->lastStep;
  server->current.frames_done = 0;
  server->current.iterations_done = 0;
  setAtoms();
  setBonds();
  //start server now initial data set
  server->startServer();
#endif
}

void OutputFAHGUI::doRun(int step) {
  GUIServer::request_t request;

  switch (request = server->getRequest()) {
  case GUIServer::GS_META_REQUEST:
  case GUIServer::GS_COORD_REQUEST:
    server->startUpdate();
          
    if (request == GUIServer::GS_META_REQUEST) {
      server->info.iterations = app->lastStep / 1000;
      server->info.frames = app->lastStep;
      setAtoms();
      setBonds();
            
    } else {
#ifdef HAVE_LIBFAH
      server->current.iterations_done = app->currentStep / 1000;
#endif
      server->current.frames_done = app->currentStep;
      server->current.energy = kineticEnergy(app->topology, &app->velocities);
      server->current.temperature =
        temperature(app->topology, &app->velocities);

      setCoords();
    }

    server->endUpdate();
    break;
          
  default: break;
  }
}

void OutputFAHGUI::doFinalize(int step) {
  doRun(step);
  if (server) delete server;
  server = 0;
}

Output *OutputFAHGUI::doMake(const vector<Value> &values) const {
  return new OutputFAHGUI(values[0], values[1], values[2], values[3], values[4]);
}

bool OutputFAHGUI::isIdDefined(const Configuration *config) const {
  return config->valid(getId());
}

void OutputFAHGUI::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(name, ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(myOutputFreq, ConstraintValueType::Positive())));
  parameter.push_back  
    (Parameter(keyword + "Port",
               Value(myPort, ConstraintValueType::Positive())));
  parameter.push_back  
    (Parameter(keyword + "PortRange",
               Value(myPortRange, ConstraintValueType::Positive())));
  parameter.push_back  
    (Parameter(keyword + "Proj",
               Value(myProjName, ConstraintValueType::NoConstraints())));
}

bool OutputFAHGUI::adjustWithDefaultParameters(vector<Value> &values,
                                               const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (!values[0].valid()) values[0] = name;
  if (!values[1].valid()) values[1] = 1;
  //if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
  //  values[1] = (*config)[InputOutputfreq::keyword];
  if (!values[2].valid()) values[2] = 52753;
  if (!values[3].valid()) values[3] = 1;
  if (!values[4].valid()) values[4] = "Protomol_3.0";

  return checkParameters(values);
}

void OutputFAHGUI::setCoords() {
  Real x, y, z, sz; 
  Vector3D posMi;

  x = y = z = 0.0;
  sz = app->positions.size();
  for (unsigned int i = 0; i < app->positions.size(); i++) {
    x += app->positions[i].c[0];
    y += app->positions[i].c[1];
    z += app->positions[i].c[2];
  }

  x /= sz; y /= sz; z /= sz;
  for (unsigned int i = 0; i < app->positions.size(); i++) {
    posMi = app->topology->minimalPosition(app->positions[i]);
    //server->xyz[i].x = app->positions[i].c[0] - x;
    //server->xyz[i].y = app->positions[i].c[1] - y;
    //server->xyz[i].z = app->positions[i].c[2] - z;
    server->xyz[i].x = posMi.c[0] - x;
    server->xyz[i].y = posMi.c[1] - y;
    server->xyz[i].z = posMi.c[2] - z;
  }
}

void OutputFAHGUI::setBonds() {
  for (unsigned int i = 0; i < app->topology->bonds.size(); i++ ) {
    // FAH requires a < b
    if (app->topology->bonds[i].atom1 < app->topology->bonds[i].atom2) {
      server->bonds[i].a = app->topology->bonds[i].atom1;
      server->bonds[i].b = app->topology->bonds[i].atom2;

    } else {
      server->bonds[i].a = app->topology->bonds[i].atom2;
      server->bonds[i].b = app->topology->bonds[i].atom1;
    }
  }
}

void OutputFAHGUI::setAtoms() {
  float radius = 0;

  for (unsigned int i = 0; i < (app->topology->atoms).size(); i++) {
    // Determine diamiter/ set name
    int atomNameLen = app->topology->atoms[i].name.length();
    for(int j=0;j<min(atomNameLen,4);j++) server->atoms[i].type[j] = (app->topology->atoms[i].name.c_str())[j];
    if(atomNameLen < 4)
        for(int j=atomNameLen;j<4;j++) server->atoms[i].type[j] = 0;
    switch(server->atoms[i].type[0]){
        case 'H':	radius = 1.2;
                    break;
        case 'C':	radius = 1.7;
                    break;
        case 'N':	radius = 1.55;
                    break;
        case 'O':	radius = 1.52;
                    break;
        case 'S':	radius = 1.85;
                    break;
        case 'P':	radius = 1.9;
                    break;
        default:	radius = 1.9;
                    break;
    }
    //add charge
    server->atoms[i].charge = app->topology->atoms[i].scaledCharge;
    server->atoms[i].radius = radius / 2.0;
  }
}

#endif // HAVE_GUI,HAVE_LIBFAH
