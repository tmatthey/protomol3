#include <protomol/type/ScalarStructure.h>

using namespace ProtoMol;
//____ ScalarStructure
ScalarStructure::ScalarStructure() :
  Proxy(), myDoVirial(true), myDoMolecularVirial(true), myDoOutput(true) {
  clear();
}

void ScalarStructure::clear() {
  for (int i = 0; i < LAST; i++)
    myTable[i] = 0.0;
}

ScalarStructure &ScalarStructure::intoAdd(const ScalarStructure &e) {
  for (int i = 0; i < LAST; i++)
    myTable[i] += e.myTable[i];

  return *this;
}

ScalarStructure &ScalarStructure::intoAssign(const ScalarStructure &e) {
  for (int i = 0; i < LAST; i++)
    myTable[i] = e.myTable[i];

  return *this;
}

ScalarStructure &ScalarStructure::intoSubtract(const ScalarStructure &e) {
  for (int i = 0; i < LAST; i++)
    myTable[i] -= e.myTable[i];

  return *this;
}

void ScalarStructure::addVirial(const Vector3D &force12,
                                const Vector3D &diff) {
  // Add virial
  Real xy = force12.x * diff.y;
  Real xz = force12.x * diff.z;
  Real yz = force12.y * diff.z;

  myTable[static_cast<int>(VIRIALXX)] += force12.x * diff.x;
  myTable[static_cast<int>(VIRIALXY)] += xy;
  myTable[static_cast<int>(VIRIALXZ)] += xz;

  myTable[static_cast<int>(VIRIALYX)] += xy;
  myTable[static_cast<int>(VIRIALYY)] += force12.y * diff.y;
  myTable[static_cast<int>(VIRIALYZ)] += yz;

  myTable[static_cast<int>(VIRIALZX)] += xz;
  myTable[static_cast<int>(VIRIALZY)] += yz;
  myTable[static_cast<int>(VIRIALZZ)] += force12.z * diff.z;
}

void ScalarStructure::addMolVirial(const Vector3D &force12,
                                   const Vector3D &diff) {
  // Add molecular virial

  Real xy = force12.x * diff.y;
  Real xz = force12.x * diff.z;
  Real yz = force12.y * diff.z;

  myTable[static_cast<int>(MOLVIRIALXX)] += force12.x * diff.x;
  myTable[static_cast<int>(MOLVIRIALXY)] += xy;
  myTable[static_cast<int>(MOLVIRIALXZ)] += xz;

  myTable[static_cast<int>(MOLVIRIALYX)] += xy;
  myTable[static_cast<int>(MOLVIRIALYY)] += force12.y * diff.y;
  myTable[static_cast<int>(MOLVIRIALYZ)] += yz;

  myTable[static_cast<int>(MOLVIRIALZX)] += xz;
  myTable[static_cast<int>(MOLVIRIALZY)] += yz;
  myTable[static_cast<int>(MOLVIRIALZZ)] += force12.z * diff.z;
}

void ScalarStructure::addVirial(const Vector3D &force12, const Vector3D &diff,
                                const Vector3D &comDiff) {
  // Add virial
  Real xy = force12.x * diff.y;
  Real xz = force12.x * diff.z;
  Real yz = force12.y * diff.z;

  myTable[static_cast<int>(VIRIALXX)] += force12.x * diff.x;
  myTable[static_cast<int>(VIRIALXY)] += xy;
  myTable[static_cast<int>(VIRIALXZ)] += xz;

  myTable[static_cast<int>(VIRIALYX)] += xy;
  myTable[static_cast<int>(VIRIALYY)] += force12.y * diff.y;
  myTable[static_cast<int>(VIRIALYZ)] += yz;

  myTable[static_cast<int>(VIRIALZX)] += xz;
  myTable[static_cast<int>(VIRIALZY)] += yz;
  myTable[static_cast<int>(VIRIALZZ)] += force12.z * diff.z;

  // Add molecular virial
  xy = force12.x * comDiff.y;
  xz = force12.x * comDiff.z;
  yz = force12.y * comDiff.z;

  myTable[static_cast<int>(MOLVIRIALXX)] += force12.x * comDiff.x;
  myTable[static_cast<int>(MOLVIRIALXY)] += xy;
  myTable[static_cast<int>(MOLVIRIALXZ)] += xz;

  myTable[static_cast<int>(MOLVIRIALYX)] += xy;
  myTable[static_cast<int>(MOLVIRIALYY)] += force12.y * comDiff.y;
  myTable[static_cast<int>(MOLVIRIALYZ)] += yz;

  myTable[static_cast<int>(MOLVIRIALZX)] += xz;
  myTable[static_cast<int>(MOLVIRIALZY)] += yz;
  myTable[static_cast<int>(MOLVIRIALZZ)] += force12.z * comDiff.z;
}

bool ScalarStructure::virial(bool doVirial) {
  bool tmp = myDoVirial;
  myDoVirial = doVirial;
  return tmp;
}

bool ScalarStructure::output(bool doOutput) {
  bool tmp = myDoOutput;
  myDoOutput = doOutput;
  return tmp;
}

bool ScalarStructure::molecularVirial(bool doMolecularVirial) {
  bool tmp = myDoMolecularVirial;
  myDoMolecularVirial = doMolecularVirial;
  return tmp;
}

bool ScalarStructure::trajectory(bool doTrajectory) {
  bool tmp = myDoTrajectory;
  myDoTrajectory = doTrajectory;
  return tmp;
}
