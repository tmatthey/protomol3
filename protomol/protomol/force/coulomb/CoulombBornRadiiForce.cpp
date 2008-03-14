#include <protomol/force/coulomb/CoulombBornRadiiForce.h>

using namespace std;
using namespace ProtoMol;
//____ CoulombForce

const string CoulombBornRadiiForce::keyword("CoulombBornRadii");
const string CoulombBornRadiiForce::C1::keyword("C1");
const string CoulombBornRadiiForce::C2::keyword("C2");
const string CoulombBornRadiiForce::C3::keyword("C3");
const string CoulombBornRadiiForce::C4::keyword("C4");

void CoulombBornRadiiForce::getParameters(vector<Parameter> &parameters) const {
}

CoulombBornRadiiForce CoulombBornRadiiForce::make(const vector<Value> &values) {
  return CoulombBornRadiiForce();
}
