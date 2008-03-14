#include "DataComparator.h"

#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>
#include <protomol/type/String.h>

#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/io/XYZReader.h>

#ifdef DATA_COMPARATOR_STANDALONE
#include <iostream>
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

Real DataComparator::compare(const string &file1, const string &file2) {
  Vector3DBlock data1;
  Vector3DBlock data2;

  read(file1, data1);
  read(file2, data2);

  return compare(data1, data2);
}

void DataComparator::read(const string &filename, Vector3DBlock &data) {
  if (!isAccessible(filename)) THROWS("Cannot access '" << filename << "'");

  string dir, base, ext;
  splitFileName(filename, dir, base, ext);

  if (ext == "dcd") {
    DCDTrajectoryReader reader(filename);

    if (!reader.tryFormat()) THROWS("Invalid DCD file '" << filename << '"');
    reader >> data;

  } else if (ext == "xyz") {
    XYZReader reader(filename);

    if (!reader.tryFormat()) THROWS("Invalid XYZ file '" << filename << '"');
    reader >> data;

  } else {
    XYZReader reader(filename);
    if (reader.tryFormat()) {
      reader >> data;

    } else {
      reader.close();
      DCDTrajectoryReader reader(filename);

      if (reader.tryFormat()) {
        reader >> data;

      } else THROWS("Unsupported file format '" << filename << "'");
    }
  }

#ifdef DATA_COMPARATOR_STANDALONE
  cout << "Vector3DBlock Size = " << data.size() << endl;
#endif // DATA_COMPARATOR_STANDALONE
}

#ifdef DATA_COMPARATOR_STANDALONE
int main(int argc, char *argv[]) {
  try {

    if (argc != 3 && argc != 4) {
      cerr << "Syntax: " << argv[0] << " <file1> <file2> [tolerance]" << endl;
      return -1;
    }
    
    Real d = DataComparator::compare(argv[1], argv[2]);
    
    if (argc == 4) {
      Real tolerance = String::parseDouble(argv[3]);
      
      if (tolerance < d)
        cout << "Files are not within tolerance" << endl;
      else cout << "Files match" << endl;
      
    } else cout << "Maximum difference: " << d << endl;

    return 0;

  } catch (const Exception &e) {
    cerr << e << endl;
  }
}

#endif // DATA_COMPARATOR_STANDALONE
