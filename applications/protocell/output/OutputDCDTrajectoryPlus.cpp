#include "OutputDCDTrajectoryPlus.h"
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include "DCDTrajectoryWriter.h"
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputDCDTrajectoryPlus
const string OutputDCDTrajectoryPlus::keyword("DCDFile");

OutputDCDTrajectoryPlus::OutputDCDTrajectoryPlus() :
  Output(), myDCD(NULL), myMinimalImage(false) {}

OutputDCDTrajectoryPlus::OutputDCDTrajectoryPlus(const string &filename, int freq,
                                         bool minimal) :
  Output(freq), myDCD(new DCDTrajectoryWriter(filename)),
  myMinimalImage(minimal) {}

OutputDCDTrajectoryPlus::~OutputDCDTrajectoryPlus() {
  if (myDCD != NULL) delete myDCD;
}

void OutputDCDTrajectoryPlus::doInitialize() {
  if (myDCD == NULL || !myDCD->open())
    THROW(string("Can not open '") + (myDCD ? myDCD->getFilename() : "") +
          "' for " + getId() + ".");

  //get data for data n header
  const int size = app->topology->bonds.size();

  std::ostringstream stm;

  //start
  stm << "<BONDS> " << size;

  //add bonds
  for( int i=0; i<size; i++ ){
      stm << " " << app->topology->bonds[i].atom1 << " "
              << app->topology->bonds[i].atom2;
  }

  //end
  stm << " </BONDS> " << size;

  //set it
  myDCD->setComment(stm.str());

}

void OutputDCDTrajectoryPlus::doRun(int) {
  const Vector3DBlock *pos =
    (myMinimalImage ? app->outputCache.minimalPositions() : &app->positions);

  if (!myDCD->write(*pos))
    THROW(string("Could not write ") + getId() + " '" +
          myDCD->getFilename() + "'.");
}

void OutputDCDTrajectoryPlus::doFinalize(int) {
  myDCD->close();
}

Output *OutputDCDTrajectoryPlus::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectoryPlus(values[0], values[1], values[2]);
}

void OutputDCDTrajectoryPlus::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(myDCD ? myDCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
}

bool OutputDCDTrajectoryPlus::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}

