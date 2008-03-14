#include <protomol/io/DCDTrajectoryReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/SystemUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____DCDTrajectoryReader

DCDTrajectoryReader::DCDTrajectoryReader() :
  Reader(ios::binary), myCoords(NULL), mySwapEndian(false) {}

DCDTrajectoryReader::DCDTrajectoryReader(const string &filename) :
  Reader(ios::binary, filename), myCoords(NULL), mySwapEndian(false) {}

DCDTrajectoryReader::~DCDTrajectoryReader() {
  if (myCoords != NULL)
    delete myCoords;
}

bool DCDTrajectoryReader::tryFormat() {
  if (!open())
    return false;

  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  file.seekg(0, ios::beg);

  int32 n = 0;
  File::read(reinterpret_cast<char *>(&n), sizeof(int32));
  char coord[5];
  File::read(coord, 4);
  coord[4] = '\0';
  close();

  int32 m = n;
  swapBytes(m);

  if (static_cast<long>(size) >= 104 &&
      (n == 84 || m == 84) && string(coord) == "CORD")
    return !file.fail();

  file.setstate(ios::failbit);
  return false;
}

bool DCDTrajectoryReader::read() {
  if (myCoords == NULL)
    myCoords = new Vector3DBlock();
  return read(*myCoords);
}

bool DCDTrajectoryReader::read(Vector3DBlock &coords) {
  if (!is_open()) {
    // First time ...
    myX.resize(0);
    myY.resize(0);
    myZ.resize(0);

    if (!open()) return false;

    file.seekg(0, ios::end);
    ios::pos_type size = file.tellg();
    file.seekg(0, ios::beg);

    int32 n = 0;
    File::read(reinterpret_cast<char *>(&n), sizeof(int32));
    char coord[5];
    File::read(coord, 4);
    coord[4] = '\0';

    int32 m = n;
    swapBytes(m);

    // Check endianess and if the header looks ok out ...
    if (static_cast<long>(size) >= 104 &&
        (n == 84 || m == 84) && string(coord) == "CORD")
      if (m == 84) {
        mySwapEndian = true;
        report << hint << "[DCDTrajectoryReader::read] Reading " <<
        (ISLITTLEENDIAN ? "big" : "little")
               << "endian input on " <<
        (ISLITTLEENDIAN ? "little" : "big") << "endian machine." << endr;
      } else {
        file.setstate(ios::failbit);
        close();
        return false;
      }

    int32 freeIndexes = 0;
    file.seekg(40, ios::beg);
    File::read(reinterpret_cast<char *>(&freeIndexes), sizeof(int32));
    if (mySwapEndian) swapBytes(freeIndexes);

    // Skip header
    file.seekg(84 + 4, ios::beg);
    n = 0;
    File::read(reinterpret_cast<char *>(&n), sizeof(int32));
    if (mySwapEndian) swapBytes(n);
    if (n != 84 || file.fail()) {
      file.setstate(ios::failbit);
      close();
      return false;
    }
    // Skip titles
    int32 l = 0;
    n = 0;
    m = 0;
    File::read(reinterpret_cast<char *>(&n), sizeof(int32));
    if (mySwapEndian) swapBytes(n);
    File::read(reinterpret_cast<char *>(&m), sizeof(int32));
    if (mySwapEndian) swapBytes(m);
    //file.seekg (m*80, ios::cur);
    comment.resize(m * 80 + m - 1);
    for (unsigned int i = 0, j = 0; i < (unsigned int)m; i++, j += 81) {
      File::read(&comment[j], 80);
      comment[j + 80] = '\n';
    }

    File::read(reinterpret_cast<char *>(&l), sizeof(int32));
    if (mySwapEndian) swapBytes(l);
    if ((n - 4) % 80 != 0 || file.fail() || l != n) {
      file.setstate(ios::failbit);
      close();
      return false;
    }

    n = 0;
    l = 0;
    int32 count = 0;
    File::read(reinterpret_cast<char *>(&n), sizeof(int32));
    if (mySwapEndian) swapBytes(n); // 4
    File::read(reinterpret_cast<char *>(&count), sizeof(int32));
    if (mySwapEndian) swapBytes(count); // number of atoms
    File::read(reinterpret_cast<char *>(&l), sizeof(int32));
    if (mySwapEndian) swapBytes(l); // 4
    if (n != 4 || l != 4 || file.fail()) {
      file.setstate(ios::failbit);
      close();
      return false;
    }

    // Skip free indexes
    if (freeIndexes > 0)
      file.seekg(4 * (count - freeIndexes + 2), ios::cur);

    myX.resize(count);
    myY.resize(count);
    myZ.resize(count);
  }

  // Read next frame
  // X-dim
  int32 count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count); // number of atoms
  count /= sizeof(int32);
  if ((unsigned int)count != myX.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  File::read(reinterpret_cast<char *>(&(myX[0])), sizeof(float4) * count);
  count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count);
  count /= sizeof(int32);
  if ((unsigned int)count != myX.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  // Y-dim
  count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count); // number of atoms
  count /= sizeof(int32);
  if ((unsigned int)count != myY.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  File::read(reinterpret_cast<char *>(&(myY[0])), sizeof(float4) * count);
  count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count);
  count /= sizeof(int32);
  if ((unsigned int)count != myY.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  // Z-dim
  count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count); // number of atoms
  count /= sizeof(int32);
  if ((unsigned int)count != myZ.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  File::read(reinterpret_cast<char *>(&(myZ[0])), sizeof(float4) * count);
  count = 0;
  File::read(reinterpret_cast<char *>(&count), sizeof(int32));
  if (mySwapEndian) swapBytes(count);
  count /= sizeof(int32);
  if ((unsigned int)count != myZ.size() || file.fail()) {
    file.setstate(ios::failbit);
    close();
    return false;
  }

  // Copy back to right structure
  coords.resize(count);
  for (int i = 0; i < count; ++i) {
    if (mySwapEndian) swapBytes(myX[i]);
    coords[i].x = myX[i];
    if (mySwapEndian) swapBytes(myY[i]);
    coords[i].y = myY[i];
    if (mySwapEndian) swapBytes(myZ[i]);
    coords[i].z = myZ[i];
  }

  return !file.fail();
}

XYZ DCDTrajectoryReader::getXYZ() const {
  XYZ res;
  if (myCoords != NULL)
    res.coords = (*myCoords);
  res.names.resize(res.coords.size(), "NONAME");
  return res;
}

Vector3DBlock *DCDTrajectoryReader::orphanCoords() {
  Vector3DBlock *tmp = myCoords;
  myCoords = NULL;
  return tmp;
}

DCDTrajectoryReader &operator>>(DCDTrajectoryReader &dcdTrajectoryReader,
                                XYZ &xyz) {
  dcdTrajectoryReader.read(xyz.coords);
  if (xyz.coords.size() != xyz.names.size())
    xyz.names.resize(xyz.coords.size(), "NONAME");
  return dcdTrajectoryReader;
}

DCDTrajectoryReader &operator>>(DCDTrajectoryReader &dcdTrajectoryReader,
                                Vector3DBlock &coords) {
  dcdTrajectoryReader.read(coords);
  return dcdTrajectoryReader;
}
