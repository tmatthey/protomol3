#include "DataComparator.h"

#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>
#include <protomol/type/String.h>

#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/io/XYZTrajectoryReader.h>
#include <protomol/io/XYZReader.h>
#include <protomol/type/XYZ.h>

#ifdef DATA_COMPARATOR_STANDALONE
#include <iostream>
#include <iomanip>
#endif

using namespace std;
using namespace ProtoMol;

Real DataComparator::compare(const Real &data1, const Real &data2) {
  return RealAbs(data1 - data2);
}

Real DataComparator::compare(const Vector3D &data1, const Vector3D &data2) {
  Real max = 0;
  Real d;

  d = compare(data1.x, data2.x);
  if (d > max) max = d;
  d = compare(data1.y, data2.y);
  if (d > max) max = d;
  d = compare(data1.z, data2.z);
  if (d > max) max = d;

  return max;
}

Real DataComparator::compare(const Vector3DBlock &data1,
                             const Vector3DBlock &data2) {
  Real max = 0;
  Real d;

  Vector3DBlock::const_iterator it1;
  Vector3DBlock::const_iterator it2;

  for (it1 = data1.begin(), it2 = data2.begin();
       it1 != data1.end() && it2 != data2.end(); it1++, it2++) {
    d = compare(*it1, *it2);
    if (d > max) max = d;
  }    

  if (it1 != data1.end() || it2 != data2.end())
    THROW("Vector3DBlock sizes don't match");

  return max;  
}

Real DataComparator::compare(const vector<XYZ> &data1,
                             const vector<XYZ> &data2, Real tolerance,
                             unsigned int &divergeFrame) {
  Real max = 0;
  Real d;

  divergeFrame = 0;

  if (data1.size() != data2.size())
    THROWS("Frame sizes not equal " << data1.size() << " != " << data2.size());

  for (unsigned int i = 0; i < data1.size(); i++) {
    d = compare(data1[i].coords, data2[i].coords);
    if (d > max) max = d;

    if (max > tolerance && !divergeFrame)
      divergeFrame = i + 1;
  }

  return max;
}

Real DataComparator::compare(const string &file1, const string &file2,
                             Real tolerance, unsigned int &divergeFrame) {
  vector<XYZ> data1;
  vector<XYZ> data2;

  read(file1, data1);
  read(file2, data2);

  return compare(data1, data2, tolerance, divergeFrame);
}

void DataComparator::read(const string &filename, vector<XYZ> &data) {
  if (!isAccessible(filename)) THROWS("Cannot access '" << filename << "'");

  XYZTrajectoryReader reader(filename);
  data.clear();

  try {
    reader >> data;

  } catch (const Exception &e1) {
    reader.close();
    XYZReader reader(filename);

    try {
      XYZ xyz;
      reader >> xyz;
      data.push_back(xyz);

    } catch (const Exception &e2) {
      reader.close();
      DCDTrajectoryReader reader(filename);

      try {
        reader >> data;
        
      } catch (const Exception &e3) {
        THROWS("Unsupported file format '" << filename << "' due to errors:"
               << endl << e1 << endl << e2 << endl << e3 << endl);
      }
    }
  }

#ifdef DATA_COMPARATOR_STANDALONE
  cout << "Frames: " << setw(5) << data.size();
  if (data.size()) cout << " Atoms: " << setw(5) << data[0].size();
  cout << endl;
#endif // DATA_COMPARATOR_STANDALONE
}

#ifdef DATA_COMPARATOR_STANDALONE
#ifndef _WIN32
#include <protomol/debug/Debugger.h>
#endif

int main(int argc, char *argv[]) {
  Exception::enableStackTraces = true;
#ifndef _WIN32
  Debugger::initStackTrace(argv[0]);
#endif

  bool result = true;

  try {

    if (argc != 3 && argc != 4) {
      cerr << "Syntax: " << argv[0] << " <file1> <file2> [tolerance]" << endl;
      return -1;
    }

    Real tolerance = 0;
    unsigned int divergeFrame = 0;

    if (argc == 4) tolerance = String::parseDouble(argv[3]);

    Real d = DataComparator::compare(argv[1], argv[2], tolerance, divergeFrame);
    
    if (argc == 4) {
      if (tolerance < d) {
        cout << "Files do not match" << endl
             << "  Tolerance       = " << tolerance << endl
             << "  Max diff        = " << d << endl
             << "  Divergent frame = " << divergeFrame << endl;

        result = false;

      } else cout << "Files match" << endl;
      
    } else {
      cout << "Maximum difference: " << d << endl;
      if (d) result = false;
    }

  } catch (const Exception &e) {
    cerr << e << endl;
    cerr << setw(10) << getFileSize(argv[1]) << " bytes " << argv[1] << endl;
    cerr << setw(10) << getFileSize(argv[2]) << " bytes " << argv[2] << endl;
  }

  return result ? 0 : 1;
}

#endif // DATA_COMPARATOR_STANDALONE
