#ifndef PROTOMOL_DATA_COMPARATOR_H
#define PROTOMOL_DATA_COMPARATOR_H

#include <string>
#include <protomol/type/Real.h>

namespace ProtoMol {
  class Vector3D;
  class Vector3DBlock;

  class DataComparator {
  public:
    static Real compare(const Real &data1, const Real &data2);
    static Real compare(const Vector3D &data1, const Vector3D &data2);
    static Real compare(const Vector3DBlock &data1,
                        const Vector3DBlock &data2);
    static Real compare(const std::string &file1, const std::string &file2);
    static void read(const std::string &filename, Vector3DBlock &data);
  };
}

#endif // PROTOMOL_DATA_COMPARATOR_H

