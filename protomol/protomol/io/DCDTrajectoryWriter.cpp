#include <protomol/io/DCDTrajectoryWriter.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____DCDTrajectoryWriter

DCDTrajectoryWriter::DCDTrajectoryWriter(Real timestep, unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc), myIsLittleEndian(isLittleEndian),
  myFirstStep(firststep), myTimeStep(timestep) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const string &filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc, filename),
  myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep) {}

DCDTrajectoryWriter::DCDTrajectoryWriter(const char *filename, Real timestep,
                                         unsigned int firststep,
                                         bool isLittleEndian) :
  Writer(ios::binary | ios::trunc, string(filename)),
  myIsLittleEndian(isLittleEndian), myFirstStep(firststep),
  myTimeStep(timestep) {}

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
  if (!is_open())
    File::open(filename.c_str(), ios::binary | ios::in | ios::out);

  // Try to read the number of frames
  file.clear();
  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  if (file.fail()) return false;

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
    file.seekg(8, ios::beg);
    read((char *)&numSets, 4);

    if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(numSets);
    ++numSets;
    if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(numSets);

    file.seekp(8, ios::beg);
    //  8: Number of sets of coordinates, NAMD=0 ???
    file.write((char *)&numSets, 4);  
    file.seekp(20, ios::beg);
    // 20: Number of sets of coordinates, NAMD=0 ???
    file.write((char *)&numSets, 4);   

    file.seekg(0, ios::end);

  } else {
    // First time ...
    close();
    open(filename.c_str(), ios::binary | ios::out | ios::trunc);

    // Write header
    file.write((char *)&n84, 4); //  0
    file.write(string("CORD").c_str(), 4); //  4
    //  8: Number of sets of coordinates, NAMD=0 ???
    file.write((char *)&numSets, 4);
    // 12: Starting timestep of DCD file, should never be zero
    file.write((char *)&firstStep, 4);
    // 16: Timesteps between DCD saves 
    file.write((char *)&numSteps, 4);
    // 20: NAMD writes += numSteps 
    file.write((char *)&numSets, 4);
    file.write((char *)&n0, 4);   // 24
    file.write((char *)&n0, 4);   // 28
    file.write((char *)&n0, 4);   // 32
    file.write((char *)&n0, 4);   // 36
    file.write((char *)&n0, 4);   // 40
    // 44 : length of a timestep
    file.write((char *)&timeStep, 4);
    // 48 : unit cell, none=0, used=1
    file.write((char *)&n0, 4);   
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    file.write((char *)&n0, 4);
    // Pretend to be Charmm 24
    file.write((char *)&n24, 4);   
    file.write((char *)&n84, 4);

    // Write DCD title record
    file.write((char *)&n164, 4);
    file.write((char *)&n2, 4);
    string remarks = string("Remarks: File '") + filename + "'. ProtoMol (" +
      __DATE__ + " at " + __TIME__ + ")";
    file.write(getRightFill(remarks, 80).c_str(), 80);
    file.write(getRightFill(string("Remarks: " + comment), 80).c_str(), 80);
    file.write((char *)&n164, 4);

    // Write DCD num-atoms record
    file.write((char *)&n4, 4);
    file.write((char *)&nAtoms, 4);
    file.write((char *)&n4, 4);
  }

  return !file.fail();
}

bool DCDTrajectoryWriter::write(const Vector3DBlock &coords) {
  const unsigned int count = coords.size();
  if (!reopen(count)) return false;

  myX.resize(count);
  myY.resize(count);
  myZ.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    myX[i] = static_cast<float>(coords[i].c[0]);
    myY[i] = static_cast<float>(coords[i].c[1]);
    myZ[i] = static_cast<float>(coords[i].c[2]);
    if (myIsLittleEndian != ISLITTLEENDIAN) {
      swapBytes(myX[i]);
      swapBytes(myY[i]);
      swapBytes(myZ[i]);
    }
  }

  int32 nAtoms = static_cast<int32>(count * 4);
  if (myIsLittleEndian != ISLITTLEENDIAN) swapBytes(nAtoms);

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myX[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myY[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

  file.write((char *)&nAtoms, sizeof(int32));
  file.write((char *)&(myZ[0]), count * sizeof(float4));
  file.write((char *)&nAtoms, sizeof(int32));

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
