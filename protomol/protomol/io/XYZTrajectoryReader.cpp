#include <protomol/io/XYZTrajectoryReader.h>

#include <protomol/type/String.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol;

//____XYZTrajectoryReader

XYZTrajectoryReader::XYZTrajectoryReader() : xyz(0) {}

XYZTrajectoryReader::XYZTrajectoryReader(const string &filename) :
  XYZReader(filename), xyz(0) {}

XYZTrajectoryReader::~XYZTrajectoryReader() {
  if (xyz) delete xyz;
}

bool XYZTrajectoryReader::tryFormat() {
  vector<XYZ> xyz;
  return read(xyz);
}

bool XYZTrajectoryReader::read() {
  if (!xyz) xyz = new vector<XYZ>();
  return read(*xyz);
}

bool XYZTrajectoryReader::read(vector<XYZ> &xyz) {
  try {
    doRead(xyz);
    return true;

  } catch (const Exception &e) {}

  return false;
}

void XYZTrajectoryReader::doRead(vector<XYZ> &xyz) {
  if (!is_open()) {
    if (!open()) THROW("open failed");

  } else file.seekg(0);

  vector<string> tokens;

  // Number of frames
  if (getLineTokens(tokens) != 1) THROW("Invalid frame count");
  unsigned int n = String::parseUInteger(tokens[0]);

  xyz.resize(n);

  try {
    // Read frames
    for (unsigned int i = 0; i < n && !file.fail(); ++i)
      XYZReader::doRead(xyz[i].coords, xyz[i].names);
    
    if (file.fail()) THROW("Reading data failed");

  } catch (const Exception &e) {
    xyz.clear();
    throw e;
  }
}

vector<XYZ> *XYZTrajectoryReader::orphanXYZ() {
  vector<XYZ> *tmp = xyz;
  xyz = 0;
  return tmp;
}

XYZTrajectoryReader &ProtoMol::operator>>(XYZTrajectoryReader &reader,
                                          vector<XYZ> &xyz) {
  reader.doRead(xyz);
  return reader;
}
