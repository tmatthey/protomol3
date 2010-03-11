#include "LangevinFlowCoupledIntegrator.h"
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ LangevinFlowCoupledIntegrator

const string LangevinFlowCoupledIntegrator::keyword("LangevinFlowCoupled");

LangevinFlowCoupledIntegrator::LangevinFlowCoupledIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1),
  averageVelocityX(0.0), averageVelocityY(0.0), averageVelocityZ(0.0)
   {}

LangevinFlowCoupledIntegrator::
LangevinFlowCoupledIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
                          int seed, Real avVX, Real avVY, Real avVZ,
                          ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed),
  averageVelocityX(avVX), averageVelocityY(avVY), averageVelocityZ(avVZ)
  {}

void LangevinFlowCoupledIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
}

void LangevinFlowCoupledIntegrator::doDrift() {
  const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinFlowCoupledIntegrator::doHalfKick() {
    const unsigned int count = app->positions.size();
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs


    Real average_velocity = 0.0;
    Real ke = 0.0;

    for (unsigned int i = 0; i < count; i++ ) {

        Vector3D fluidVelocity(averageVelocityX, averageVelocityY, averageVelocityZ);
        //Vector3D projectedVelocity(0.0, 0.0, 0.0);

        /*const unsigned int myCellCenter = app->topology->atoms[i].cell_center;
        const Real normFluidV = fluidVelocity.norm();
        Real factor = 1.0;

        //dont use velocity if center
        if(myCellCenter != i && normFluidV != 0.0){
          Vector3D sc =  app->positions[myCellCenter] - app->positions[i];
                        //app->topology->minimalDifference( app->positions[i], app->positions[myCellCenter]);
          factor = sc.dot(fluidVelocity) / ( normFluidV * sc.norm() );
          //if(factor > 1.0) report << hint << "factor too big " << factor << endr;
          if(factor > 0.0){
            //report << hint << "factor " << factor << endr;
            projectedVelocity = fluidVelocity;// * fabs(factor);
          }else{
            factor = 0.1;
          }
        }
        if( myCellCenter == i ) factor = 0.1;*/

        Real aGamma = myGamma;// * factor;  //####Removed factor
        if(aGamma < 0.1){ //####was 0.1
          aGamma = 1.0;
          //report << hint << "aGamma < 0" << endr;
        }
        const Real fdt = ( 1.0 - exp( -0.5 * aGamma * dt ) ) / aGamma;
        const Real vdt = exp(-0.5*aGamma*dt);
        const Real ndt = sqrt( ( 1.0 - exp( -aGamma * dt ) ) / (2.0 * aGamma) );
        const Real forceConstant = 2 * Constant::BOLTZMANN * myLangevinTemperature *
                  aGamma; //SI::BOLTZMANN* 1.0e21

        //  Generate gaussian random numbers for each spatial direction
        Vector3D gaussRandCoord1(randomGaussianNumber(mySeed),
                                 randomGaussianNumber(mySeed),
                                 randomGaussianNumber(mySeed));
        Real mass = app->topology->atoms[i].scaledMass;
        Real sqrtFCoverM = sqrt(forceConstant / mass);
        //Real sqrtFCoverMf = sqrt(forceConstant) / mass;
        // fluid velocity (use fixed initially with random purtubation)
        //Vector3D fluidVelocity(0.0 + (randomNumber(1234) - 0.5) * 0.1,0.00+(randomNumber(1234) - 0.5) * 0.1,(randomNumber(1234) - 0.5) * 0.1),
        //                       projectedVelocity(0,0,0);

        //remove velocity
        app->velocities[i] -= fluidVelocity;

        // semi-update velocities
        app->velocities[i] = app->velocities[i]*vdt
                                +(*myForces)[i] * fdt / mass
                                  //+ projectedVelocity * fdt
                                    +gaussRandCoord1*sqrtFCoverM*ndt;
                                    //+gaussRandCoord1*(sqrt(sqrtFCoverM + sqrtFCoverMf))*ndt;

        //find "real" temperature
        for(int k=0; k<3; k++){
          ke += 0.5 * app->velocities[i].c[k] * app->velocities[i].c[k] * mass;
        }

        //replace velocity
        app->velocities[i] += fluidVelocity;

        //find average
        average_velocity += app->velocities[i].c[0];



    }
    //report << hint << "Average velocity  " << average_velocity / (Real)count
    //        << " Set velocity " << averageVelocityX << " Temp " << 2.0 * ke / Constant::BOLTZMANN / count / 3.0 << endr;

    buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinFlowCoupledIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature", Value(myLangevinTemperature,
                                    ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234));
  parameters.push_back
    (Parameter("averageflowvelocityx", Value(averageVelocityX, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity x")));
  parameters.push_back
    (Parameter("averageflowvelocityy", Value(averageVelocityY, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity y")));
  parameters.push_back
    (Parameter("averageflowvelocityz", Value(averageVelocityZ, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity z")));
}

STSIntegrator *LangevinFlowCoupledIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinFlowCoupledIntegrator(values[0], values[1], values[2],
                                            values[3], values[4], values[5], values[6],
                                            fg);
}
