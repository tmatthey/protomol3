#include <protomol/io/XYZTrajectoryReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____XYZTrajectoryReader

XYZTrajectoryReader::XYZTrajectoryReader() :
  Reader(), myCoords(NULL), myNames(NULL) {}

XYZTrajectoryReader::XYZTrajectoryReader(const string &filename) :
  Reader(filename), myCoords(NULL), myNames(NULL) {}

XYZTrajectoryReader::~XYZTrajectoryReader() {
  if (myCoords != NULL) delete myCoords;
  if (myNames != NULL) delete myNames;
}

bool XYZTrajectoryReader::tryFormat() {
  if (!open())
    return false;

  // Number of frames and atoms
  unsigned int frames, count = 0;
  file >> frames >> count;

  // Atoms
  string str;
  Real x;
  file >> str >> x >> x >> x;
  close();
  return !file.fail();
}

bool XYZTrajectoryReader::read() {
  if (myCoords == NULL)
    myCoords = new Vector3DBlock();
  if (myNames == NULL)
    myNames = new vector<string>();
  return read(*myCoords, *myNames);
}

bool XYZTrajectoryReader::read(XYZ &xyz) {
  return read(xyz.coords, xyz.names);
}

bool XYZTrajectoryReader::read(Vector3DBlock &coords, vector<string> &names) {
  if (!is_open()) {
    if (!open()) return false;

    // Number of frames;
    unsigned int n = 0;
    file >> n;
    if (n == 0 && file.good())
      return true;
  }

  // Number of atoms;
  unsigned int n = 0;
  file >> n;
  if (file.fail()) {
    close();
    return false;
  }

  coords.resize(n);
  names.resize(n);

  // Read atoms
  for (unsigned int i = 0; i < n && !file.fail(); ++i)
    file >> names[i] >> coords[i].x >> coords[i].y >> coords[i].z;

  return !file.fail();
}

XYZ XYZTrajectoryReader::getXYZ() const {
  XYZ res;
  if (myCoords != NULL)
    res.coords = (*myCoords);
  if (myNames != NULL)
    res.names = (*myNames);
  return res;
}

Vector3DBlock *XYZTrajectoryReader::orphanCoords() {
  Vector3DBlock *tmp = myCoords;
  myCoords = NULL;
  return tmp;
}

vector<string> *XYZTrajectoryReader::orphanNames() {
  vector<string> *tmp = myNames;
  myNames = NULL;
  return tmp;
}

XYZTrajectoryReader &ProtoMol::operator>>(XYZTrajectoryReader &xyzReader,
                                          XYZ &xyz) {
  xyzReader.read(xyz.coords, xyz.names);
  return xyzReader;
}

XYZTrajectoryReader &ProtoMol::operator>>(XYZTrajectoryReader &xyzReader,
                                          Vector3DBlock &coords) {
  if (xyzReader.myNames == NULL)
    xyzReader.myNames = new vector<string>();
  xyzReader.read(coords, *xyzReader.myNames);
  return xyzReader;
}
