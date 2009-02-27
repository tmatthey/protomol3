#include <protomol/integrator/base/OpenMMIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Zap.h>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ OpenMMIntegrator

const string OpenMMIntegrator::keyword("OpenMM");

OpenMMIntegrator::OpenMMIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1),
  HarmonicBondForce(false), HarmonicAngleForce(false), NonbondedForce(false)
  {
#if defined (HAVE_OPENMM)
    system = 0;
    bonds = 0;
    angles = 0;    
    nonbonded = 0;
    integrator = 0;
    context = 0;
#endif
  }

OpenMMIntegrator::
OpenMMIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
                          int seed,
                          bool bond, bool angle, bool nonbond,
                          ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed),
  HarmonicBondForce(bond), HarmonicAngleForce(angle), NonbondedForce(nonbond)
  {
#if defined (HAVE_OPENMM)
    system = 0;
    bonds = 0;
    angles = 0;
    nonbonded = 0;
    integrator = 0;
    context = 0;
#endif
  }

OpenMMIntegrator::~OpenMMIntegrator() {
#if defined (HAVE_OPENMM)
  zap(context);
  zap(integrator);
  //zap(nonbonded);
  //zap(angles);
  zap(system);
  //zap(bonds);
#endif
}

void OpenMMIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  //openMM

#if defined (HAVE_OPENMM)
  unsigned int sz = app->positions.size();


  system = new OpenMM::System(sz, 0);
  for (int i = 0; i < sz; ++i)
    system->setParticleMass(i, app->topology->atoms[i].scaledMass);

  //openMM forces

  if ( HarmonicBondForce ){
    unsigned int numBonds = app->topology->bonds.size();
    bonds = new OpenMM::HarmonicBondForce(numBonds);
    system->addForce(bonds);

    for (int i = 0; i < numBonds; ++i){
      unsigned int a1 = app->topology->bonds[i].atom1; unsigned int a2 = app->topology->bonds[i].atom2;
      Real r_0 = app->topology->bonds[i].restLength  * Constant::ANGSTROM_NM;
      Real k = app->topology->bonds[i].springConstant  * Constant::KCAL_KJ * 2.0; //times 2 as Amber is 1/2 k(b-b_0)^2
      bonds->setBondParameters(i, a1, a2, r_0, k);
    }
  }

  if ( HarmonicAngleForce ){
    unsigned int numAngles = app->topology->angles.size();
    angles = new OpenMM::HarmonicAngleForce(numAngles);
    system->addForce(angles);
    for (int i = 0; i < numAngles; i++){
        unsigned int a1 = app->topology->angles[i].atom1;
        unsigned int a2 = app->topology->angles[i].atom2;
        unsigned int a3 = app->topology->angles[i].atom3;
        Real theta0 = acos(cos(app->topology->angles[i].restAngle));
        Real k_t = app->topology->angles[i].forceConstant * Constant::KCAL_KJ / 50.0;
        report << hint << "rest angle " << theta0 << " k " << k_t << " ang size " << numAngles << endr;
        angles->setAngleParameters(i, a3, a2, a1, theta0, k_t);
    }
  }

  if ( NonbondedForce ){
    nonbonded = new OpenMM::NonbondedForce(sz,0);
    system->addForce(nonbonded);
    for (int i = 0; i < sz; i++){
      int type1 = app->topology->atoms[i].type;
      //report << hint << "sigma " <<  app->topology->atomTypes[type1].sigma << ", eps " << app->topology->atomTypes[type1].epsilon << endr;
      //nonbonded->setParticleParameters(i, 0, //app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT, 
          //app->topology->atomTypes[type1].sigma * Constant::ANGSTROM_NM, app->topology->atomTypes[type1].epsilon * -Constant::KCAL_KJ);
    }
  }

  //openMM Initialize
  integrator = new OpenMM::LangevinIntegrator(myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS); //ps
  context = new OpenMM::OpenMMContext(*system, *integrator);

  OpenMM::Vec3 openMMvecp, openMMvecv;
  for (int i = 0; i < sz; ++i){
    for (int j = 0; j < 3; j++){
      openMMvecp[j] = app->positions[i].c[j] * Constant::ANGSTROM_NM;
      openMMvecv[j] = app->velocities[i].c[j] * Constant::ANGSTROM_NM * Constant::INV_TIMEFACTOR;
    }
    openMMpositions.push_back(openMMvecp);
    openMMvelocities.push_back(openMMvecv);
  }

  context->setPositions(openMMpositions);
  context->setVelocities(openMMvelocities);

  //print platform
  report << plain << "OpenMM platform is: '" << context->getPlatform().getName() << "'." << endr;

#else

  //print platform
  report << plain << "OpenMM platform is not available." << endr;

#endif

}

void OpenMMIntegrator::run(int numTimesteps) {

  preStepModify();

#if defined (HAVE_OPENMM)
  unsigned int sz = app->positions.size();

  // do integration
  integrator->step(numTimesteps);

  // Retrive data
  const OpenMM::State state = context->getState(OpenMM::State::Positions | 
                                                OpenMM::State::Velocities |
                                                OpenMM::State::Forces |
                                                OpenMM::State::Energy);
  openMMpositions = state.getPositions();
  openMMvelocities = state.getVelocities();
  openMMforces = state.getForces();

  for (int i = 0; i < sz; ++i){
   for (int j = 0; j < 3; j++){
     app->positions[i].c[j] = openMMpositions[i][j] * Constant::NM_ANGSTROM; //nm to A
     app->velocities[i].c[j] = openMMvelocities[i][j];// * Constant::NM_ANGSTROM * Constant::TIMEFACTOR; //nm/ps to A/fs?
     (*myForces)[i].c[j] = openMMforces[i][j] * Constant::NM_ANGSTROM * Constant::KJ_KCAL; //KJ/nm to Kcal/A
    }
  }

#endif

  postStepModify();

}

void OpenMMIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("temperature", Value(myLangevinTemperature, ConstraintValueType::NotNegative())));
  parameters.push_back(Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR), ConstraintValueType::NotNegative())));
  parameters.push_back(Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()), 1234));
  //OpenMM forces
  parameters.push_back(Parameter( "HarmonicBondForce", Value( HarmonicBondForce, ConstraintValueType::NoConstraints() ), false ));
  parameters.push_back(Parameter( "HarmonicAngleForce", Value( HarmonicAngleForce, ConstraintValueType::NoConstraints() ), false ));
  parameters.push_back(Parameter( "NonbondedForce", Value( NonbondedForce, ConstraintValueType::NoConstraints() ), false ));

}

STSIntegrator *OpenMMIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new OpenMMIntegrator(values[0], values[1], values[2],
                                       values[3], 
                                       values[4], values[5], values[6],
                                       fg);
}

// Constants for GROMACS-PROTOMOL conversion
namespace ProtoMol {
  namespace Constant {

    const Real ANGSTROM_NM = 0.1;
    const Real NM_ANGSTROM = 10.0;
    const Real KCAL_KJ = 4.184;
    const Real KJ_KCAL = 1.0 / 4.184;
    const Real FS_PS = 0.001;

  }
}
