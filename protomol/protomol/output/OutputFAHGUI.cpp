#ifdef HAVE_LIBFAH
#include <protomol/output/OutputFAHGUI.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#include <fah/core/GUIServer.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace FAH;

//____ OutputFAHGUI
const string OutputFAHGUI::keyword("FAHGUI");

OutputFAHGUI::OutputFAHGUI() : name("ProtoMol"), server() {}

OutputFAHGUI::OutputFAHGUI(const string &name, int freq) :
  Output(freq), name(name), server(0) {}

void OutputFAHGUI::doInitialize() {
  server = new GUIServer(name.c_str(), app->topology->atoms.size(),
                         app->topology->bonds.size());
}

void OutputFAHGUI::doRun(int step) {
  GUIServer::request_t request;

  switch (request = server->getRequest()) {
  case GUIServer::GS_META_REQUEST:
  case GUIServer::GS_COORD_REQUEST:
    server->startUpdate();
          
    if (request == GUIServer::GS_META_REQUEST) {
      server->info.iterations = app->lastStep / 1000;
      setAtoms();
      setBonds();
            
    } else {
      server->current.iterations_done = app->currentStep / 1000;
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
  return new OutputFAHGUI(values[0], values[1]);
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
}

bool OutputFAHGUI::adjustWithDefaultParameters(vector<Value> &values,
                                               const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (!values[0].valid()) values[0] = name;
  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  return checkParameters(values);
}

void OutputFAHGUI::setCoords() {
  Real x, y, z, sz; 

  x = y = z = 0.0;
  sz = app->positions.size();
  for (unsigned int i = 0; i < app->positions.size(); i++) {
    x += app->positions[i].x;
    y += app->positions[i].y;
    z += app->positions[i].z;
  }

  x /= sz; y /= sz; z /= sz;
  for (unsigned int i = 0; i < app->positions.size(); i++) {
    server->xyz[i].x = app->positions[i].x - x;
    server->xyz[i].y = app->positions[i].y - y;
    server->xyz[i].z = app->positions[i].z - z;
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
  int element;
  float mass;
  float radius = 0;

  for (unsigned int i = 0; i < (app->topology->atoms).size(); i++) {
    // Determine element by mass
    // This will not work with united atom models or with heavy hydrogens etc
    element = 0; // unknown
    mass = app->topology->atoms[i].scaledMass;

    if (mass < 1.2 && mass >= 1.0) { // hydrogen
      element = 1;
      server->atoms[i].type[0] = 'H';
      radius = 1.2;

    } else if (mass > 11.8 && mass < 12.2) { // carbon
      element = 6;
      server->atoms[i].type[0] = 'C';
      radius = 1.7;

    } else if (mass > 14.0 && mass < 15) { // nitrogen
      element = 7;
      server->atoms[i].type[0] = 'N';
      radius = 1.55;

    } else if (mass > 15.5 && mass < 16.5) { // oxygen
      element = 8;
      server->atoms[i].type[0] = 'O';
      radius = 1.52;

    } else if (mass > 31.5 && mass < 32.5) { // sulphur
      element = 16;
      server->atoms[i].type[0] = 'S';
      radius = 1.85;

    } else if (mass > 29.5 && mass < 30.5) { // phosphorus
      element = 15;
      server->atoms[i].type[0] = 'P';
      radius = 1.9;
    }

    server->atoms[i].charge = app->topology->atoms[i].scaledCharge;
    server->atoms[i].radius = radius / 2.0;
  }
}

#endif // HAVE_LIBFAH
