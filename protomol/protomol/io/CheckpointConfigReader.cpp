#include <protomol/io/CheckpointConfigReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


CheckpointConfigReader::CheckpointConfigReader() : Reader() {}


CheckpointConfigReader::CheckpointConfigReader(const string &filename) :
  Reader(filename) {}


bool CheckpointConfigReader::tryFormat() {
  if (!open()) return false;

  string header = getline();
    
  cout << header << endl;

  return header == "!Checkpoint File!";
}


bool CheckpointConfigReader::readBase(Configuration &conf, Random &rand) {
  if (!tryFormat()) {
    cout << "Invalid checkpoint" << endl;
    return false;
  }
    
  cout << "Reading checkpoint" << endl;

  int id = 0, step = 0;
    
  string line;
  while (std::getline(file, line)) {
    if (line.find("#ID") != string::npos) file >> id;
    if (line.find("#Step") != string::npos) file >> step;
    if (line.find("#Random") != string::npos) file >> rand;
  }

  // Update initial checkpoint perameters
  conf["CheckpointStart"] = id + 1;

  // Update position file
  conf["posfile"] = Append(conf["CheckpointPosBase"], id) + ".pos";

  // Update velocities file
  conf["velfile"] = Append(conf["CheckpointVelBase"], id) + ".vel";

  // Update energy file
  if (conf.valid("allEnergiesFile"))
    conf["allEnergiesFile"] = Append(conf["allEnergiesFile"], id);

  // Update firststep
  unsigned firststep = toInt(conf["firststep"]);
  conf["firststep"] = firststep + step;

  // Update total steps
  conf["numsteps"] = toInt(conf["numsteps"]) - firststep;

  return !file.fail();
}


bool CheckpointConfigReader::readIntegrator(Integrator *integ) {
  if (!tryFormat()) return false;
    
  string line;
  while (std::getline(file, line))
    if (line.find("#Integrator") != string::npos)
      file >> *integ;

  return !file.fail();
}
