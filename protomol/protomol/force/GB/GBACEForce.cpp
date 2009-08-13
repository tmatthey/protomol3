#include <protomol/force/GB/GBACEForce.h>

using namespace std;
using namespace ProtoMol;

const string GBACEForce::keyword("GBACEForce");

void GBACEForce::getParameters(vector<Parameter> &parameters) const {

     parameters.push_back(
        Parameter("-solvationparam",Value(sigma, ConstraintValueType::NoConstraints()),3.0, Text("solvation parameter")));

     parameters.push_back(
        Parameter("-watersphereradius",Value(rho_s, ConstraintValueType::NoConstraints()), 3.0, Text("solvation parameter")));

}

GBACEForce GBACEForce::make(const vector<Value> &values)
{
   return (GBACEForce(values[0], values[1]));
}
