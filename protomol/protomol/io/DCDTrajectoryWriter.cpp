#include <protomol/io/DCDTrajectoryWriter.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____DCDTrajectoryWriter

DCDTrajectoryWriter::DCDTrajectoryWriter(Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc), myIsLittleEndian(isLittleEndian),
  myFirstStep(firststep), myTimeStep(timestep) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const string &filename,
                                         Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc,
         filename), myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const char *filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc,
         string(filename)), myIsLittleEndian(isLittleEndian),
  myFirstStep(firststep), myTimeStep(timestep) {}

bool DCDTrajectoryWriter::openWith(Real timestep, unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open();
}

bool DCDTrajectoryWriter::openWith(const string &filename, Real timestep,
                                   unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open(filename);
}

bool DCDTrajectoryWriter::openWith(const char *filename, Real timestep,
                                   unsigned int firststep,
                                   bool isLittleEndian) {
  setTimestep(timestep);
  setFirststep(firststep);
  setLittleEndian(isLittleEndian);
  return open(string(filename));
}

bool DCDTrajectoryWriter::reopen(unsigned int numAtoms) {
  if (is_open()) close();

  // Try to read the number of frames
  file.clear();
  File::open(filename.c_str(), ios::binary | ios::in);
  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  close();
  if (file.fail())
    return false;

  int32 nAtoms = static_cast<int32>(numAtoms);
  int32 numSets = 1;
  int32 numSteps = 1;
  int32 firstStep = static_cast<int32>(myFirstStep);
  float4 timeStep = static_cast<float4>(myTimeStep) *
                    Constant::INV_TIMEFACTOR;

  int32 n0 = 0;
  int32 n2 = 2;
  int32 n4 = 4;
  int32 n24 = 24;
  int32 n84 = 84;
  int32 n164 = 164;

  if (myIsLittleEndian != ISLITTLEENDIAN) {
    swapBytes(nAtoms);
    swapBytes(numSets);
    swapBytes(numSteps);
    swapBytes(firstStep);
    swapBytes(timeStep);

    swapBytes(n0);
    swapBytes(n2);
    swapBytes(n4);
    swapBytes(n24);
    swapBytes(n84);
    swapBytes(n164);
  }

  if (size > static_cast<ios::pos_type>(100)) {
    // Ok, we have already written frames
    file.clear();
    open(filename.c_str(), ios::binary | ios::in);
    file.seekg(8, ios::beg);
    read(reinterpret_cast<char *>(&numSets), 4);
    close();

    if (myIsLittleEndian != ISLITTLEENDIAN)
      swapBytes(numSets);
    ++numSets;
    if (myIsLittleEndian != ISLITTLEENDIAN)
      swapBytes(numSets);

    file.clear();
    open(
      filename.c_str(), ios::binary | ios::in | ios::out);
    file.seekp(8, ios::beg);
    //  8: Number of sets of coordinates, NAMD=0 ???
    file.write(reinterpret_cast<char *>(&numSets), 4);  
    file.seekp(20, ios::beg);
    // 20: Number of sets of coordinates, NAMD=0 ???
    file.write(reinterpret_cast<char *>(&numSets), 4);   
    close();
  } else {
    // First time ...
    file.clear();
    open(
      filename.c_str(), ios::binary | ios::out | ios::trunc);

    file.write(reinterpret_cast<char *>(&n84), 4);   //  0
    file.write(string("CORD").c_str(), 4);              //  4
    //  8: Number of sets of coordinates, NAMD=0 ???
    file.write(reinterpret_cast<char *>(&numSets), 4);   
    // 12: Starting timestep of DCD file, should be never zero
    file.write(reinterpret_cast<char *>(&firstStep), 4);  
    // 16: Timesteps between DCD saves 
    file.write(reinterpret_cast<char *>(&numSteps), 4);  
    // 20: NAMD writes += numSteps 
    file.write(reinterpret_cast<char *>(&numSets), 4);   
    file.write(reinterpret_cast<char *>(&n0), 4);   // 24
    file.write(reinterpret_cast<char *>(&n0), 4);   // 28
    file.write(reinterpret_cast<char *>(&n0), 4);   // 32
    file.write(reinterpret_cast<char *>(&n0), 4);   // 36
    file.write(reinterpret_cast<char *>(&n0), 4);   // 40
    // 44 : length of a timestep
    file.write(reinterpret_cast<char *>(&timeStep), 4);   
    // 48 : unit cell, none=0, used=1
    file.write(reinterpret_cast<char *>(&n0), 4);   

    file.write(reinterpret_cast<char *>(&n0), 4);
    file.write(reinterpret_cast<char *>(&n0), 4);
    file.write(reinterpret_cast<char *>(&n0), 4);

    file.write(reinterpret_cast<char *>(&n0), 4);
    file.write(reinterpret_cast<char *>(&n0), 4);
    file.write(reinterpret_cast<char *>(&n0), 4);

    file.write(reinterpret_cast<char *>(&n0), 4);
    file.write(reinterpret_cast<char *>(&n0), 4);
    // Pretend to be Charmm 24
    file.write(reinterpret_cast<char *>(&n24), 4);   

    file.write(reinterpret_cast<char *>(&n84), 4);

    // Write DCD title record
    file.write(reinterpret_cast<char *>(&n164), 4);
    file.write(reinterpret_cast<char *>(&n2), 4);
    file.write(getRightFill(string("Remarks: File \'" +
                                   filename + "\'. ProtoMol ("
                                   + string(__DATE__) + " at "
                                   + string(__TIME__) + ")")
                            , 80).c_str(), 80);
    file.write(getRightFill(string("Remarks: " + comment),
                            80).c_str(), 80);
    file.write(reinterpret_cast<char *>(&n164), 4);

    // Write DCD num-atoms record
    file.write(reinterpret_cast<char *>(&n4), 4);
    file.write(reinterpret_cast<char *>(&nAtoms), 4);
    file.write(reinterpret_cast<char *>(&n4), 4);
    if (file.fail()) {
      close();
      return false;
    }
    close();
  }

  file.clear();
  open(
    filename.c_str(), ios::binary | ios::out | ios::app);
  return !file.fail();
}

bool DCDTrajectoryWriter::write(const Vector3DBlock &coords) {
  const unsigned int count = coords.size();
  if (!reopen(count))
    return false;

  myX.resize(count);
  myY.resize(count);
  myZ.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    myX[i] = static_cast<float>(coords[i].x);
    myY[i] = static_cast<float>(coords[i].y);
    myZ[i] = static_cast<float>(coords[i].z);
    if (myIsLittleEndian != ISLITTLEENDIAN) {
      swapBytes(myX[i]);
      swapBytes(myY[i]);
      swapBytes(myZ[i]);
    }
  }

  int32 nAtoms = static_cast<int32>(count * 4);
  if (myIsLittleEndian != ISLITTLEENDIAN)
    swapBytes(nAtoms);
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));
  file.write(reinterpret_cast<char *>(&(myX[0])), count * sizeof(float4));
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));
  file.write(reinterpret_cast<char *>(&(myY[0])), count * sizeof(float4));
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));
  file.write(reinterpret_cast<char *>(&(myZ[0])), count * sizeof(float4));
  file.write(reinterpret_cast<char *>(&nAtoms), sizeof(int32));

  close();
  return !file.fail();
}

void DCDTrajectoryWriter::setLittleEndian(bool littleEndian) {
  myIsLittleEndian = littleEndian;
}

void DCDTrajectoryWriter::setTimestep(Real timestep) {
  myTimeStep = timestep;
}

void DCDTrajectoryWriter::setFirststep(unsigned int firststep) {
  myFirstStep = firststep;
}

DCDTrajectoryWriter &ProtoMol::operator<<(DCDTrajectoryWriter &dcdWriter,
                                          const Vector3DBlock &coords) {
  dcdWriter.write(coords);
  return dcdWriter;
}

DCDTrajectoryWriter &ProtoMol::operator<<(DCDTrajectoryWriter &dcdWriter,
                                          const XYZ &xyz) {
  dcdWriter.write(xyz.coords);
  return dcdWriter;
}
