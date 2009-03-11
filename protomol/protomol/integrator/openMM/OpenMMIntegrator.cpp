#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Zap.h>
#include <protomol/topology/LennardJonesParameters.h>

#include <vector>
#include <algorithm>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ OpenMMIntegrator

const string OpenMMIntegrator::keyword("OpenMM");

const Real OpenMMIntegrator::FudgeQQ = 0.8333;
const Real OpenMMIntegrator::FudgeLJ = 0.5;

OpenMMIntegrator::OpenMMIntegrator() :
  STSIntegrator()
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
OpenMMIntegrator(Real timestep, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces)
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

struct NBForce{
	int atom1;
	int atom2;
	Real charge;
	Real sigma;
	Real epsilon;
	Real c6;
	Real c12;
	
	NBForce( int a, int b, Real c, Real s, Real e, Real cs, Real ct ){
		atom1 = a;
		atom2 = b;
		charge = c;
		sigma = s;
		epsilon = e;
		c6 = cs;
		c12 = ct;
	}
	
	bool operator< ( const NBForce &other ) const{
		if ( atom1 < other.atom1 ){
			return true;
		}
		
		if ( atom1 == other.atom1 ){
			if ( atom2 < other.atom2 ){
				return true;
			}
		}
		
		return false;
	}
};

void OpenMMIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  //openMM

#if defined (HAVE_OPENMM)
  unsigned int sz = app->positions.size();

#ifdef DEBUG
  std::ofstream mFile ( "output.txt" );
#endif

  system = new OpenMM::System(sz, 0);
  for (unsigned int i = 0; i < sz; ++i)
    system->setParticleMass(i, app->topology->atoms[i].scaledMass);

  //openMM forces
#ifdef DEBUG
  mFile << "Bonds" << std::endl;
#endif
  if ( HarmonicBondForce ){
    unsigned int numBonds = app->topology->bonds.size();
    bonds = new OpenMM::HarmonicBondForce(numBonds);
    system->addForce(bonds);

    for (unsigned int i = 0; i < numBonds; ++i){
		unsigned int a1 = app->topology->bonds[i].atom1; unsigned int a2 = app->topology->bonds[i].atom2;
		Real r_0 = app->topology->bonds[i].restLength  * Constant::ANGSTROM_NM;
		Real k = app->topology->bonds[i].springConstant  
              * Constant::KCAL_KJ * Constant::INV_ANGSTROM_NM * Constant::INV_ANGSTROM_NM * 2.0; //times 2 as Amber is 1/2 k(b-b_0)^2
	  
#ifdef DEBUG
	  if( (app->topology->atoms[ app->topology->bonds[i].atom1 ].name[0] != 'H') &&
		  (app->topology->atoms[ app->topology->bonds[i].atom2 ].name[0] != 'H') ){
		  mFile << a1 << " " << a2 << " " << r_0 << " " << k << std::endl;
	  }
#endif

	  bonds->setBondParameters(i, a1, a2, r_0, k);
    }
  }

#ifdef DEBUG
  mFile << std::endl;

  mFile << "Angles" << std::endl;
#endif

  if ( HarmonicAngleForce ){
    unsigned int numAngles = app->topology->angles.size();
    angles = new OpenMM::HarmonicAngleForce(numAngles);
    system->addForce(angles);
    for (unsigned int i = 0; i < numAngles; i++){
        unsigned int a1 = app->topology->angles[i].atom1;
        unsigned int a2 = app->topology->angles[i].atom2;
        unsigned int a3 = app->topology->angles[i].atom3;
        Real theta0 = acos(cos(app->topology->angles[i].restAngle));
        Real k_t = app->topology->angles[i].forceConstant * Constant::KCAL_KJ * 2.0; //times 2 as Amber is 1/2 k(a-a_0)^2
        //report << hint << "rest angle " << theta0 << " k " << k_t << " ang size " << numAngles << endr;

#ifdef DEBUG
        mFile << a1 << " " << a2 << " " << a3 << " " << theta0 << " " << k_t << std::endl;        
#endif

        angles->setAngleParameters(i, a1, a2, a3, theta0, k_t);
    }
  }

#ifdef DEBUG
  mFile << std::endl;

  mFile << "Periodic Force" << std::endl;
#endif

  if ( PeriodicTorsion ){

    unsigned int numPTor = app->topology->dihedrals.size();
    unsigned int totalNumPTor = 0;
    for (unsigned int i = 0; i < numPTor; i++) 
      totalNumPTor += app->topology->dihedrals[i].multiplicity;

    PTorsion = new OpenMM::PeriodicTorsionForce(totalNumPTor);//numPTor);//
    system->addForce(PTorsion);
    for (unsigned int i = 0; i < numPTor; i++){
        unsigned int a1 = app->topology->dihedrals[i].atom1;
        unsigned int a2 = app->topology->dihedrals[i].atom2;
        unsigned int a3 = app->topology->dihedrals[i].atom3;
        unsigned int a4 = app->topology->dihedrals[i].atom4;

        unsigned int multiplicity = app->topology->dihedrals[i].multiplicity;

        for (unsigned int j = 0; j < multiplicity; j++){

          unsigned int mult = app->topology->dihedrals[i].periodicity[j];
          Real phiA = app->topology->dihedrals[i].phaseShift[j];
          Real cpA = app->topology->dihedrals[i].forceConstant[j] * Constant::KCAL_KJ;

          //idef.iparams[type].pdihs.mult, idef.iparams[type].pdihs.phiA*M_PI/180.0, idef.iparams[type].pdihs.cpA
#ifdef DEBUG
          mFile << a1 << " " << a2 << " " << a3 << " " << a4 << " " 
              << mult << " " << phiA << " " << cpA << std::endl;   
#endif

          PTorsion->setTorsionParameters(i, a1, a2, a3, a4, mult, phiA, cpA);
        }
    }
  }

#ifdef DEBUG
  mFile << std::endl;

  mFile << "RBDihedrals" << std::endl;
#endif

  if ( RBDihedralForce ){
    unsigned int numRBDih = app->topology->rb_dihedrals.size();
    RBDihedral = new OpenMM::RBTorsionForce(numRBDih);
    system->addForce(RBDihedral);
    for (unsigned int i = 0; i < numRBDih; i++){
        unsigned int a1 = app->topology->rb_dihedrals[i].atom1;
        unsigned int a2 = app->topology->rb_dihedrals[i].atom2;
        unsigned int a3 = app->topology->rb_dihedrals[i].atom3;
        unsigned int a4 = app->topology->rb_dihedrals[i].atom4;
        Real C0 = app->topology->rb_dihedrals[i].C0 * Constant::KCAL_KJ;
        Real C1 = app->topology->rb_dihedrals[i].C1 * Constant::KCAL_KJ;
        Real C2 = app->topology->rb_dihedrals[i].C2 * Constant::KCAL_KJ;
        Real C3 = app->topology->rb_dihedrals[i].C3 * Constant::KCAL_KJ;
        Real C4 = app->topology->rb_dihedrals[i].C4 * Constant::KCAL_KJ;
        Real C5 = app->topology->rb_dihedrals[i].C5 * Constant::KCAL_KJ;
        if(C0 != 0 || C1 != 0 ||C2 != 0 ||C3 != 0 ||C4 != 0 ||C5 != 0 ) {

#ifdef DEBUG
          mFile << a1 << " " << a2 << " " << a3 << " " << a4 << " " 
            << C0 << " " << C1 << " " << C2 << " " << C3 << " " << C4 << " " << C5 << std::endl;     
#endif

        }
        RBDihedral->setTorsionParameters(i, a1, a2, a3, a4, C0, C1, C2, C3, C4, C5);
    }
  }

#ifdef DEBUG
  mFile << std::endl;

  mFile << "Lenards Jones Force" << std::endl;
#endif

  if ( NonbondedForce ){

    //get 1-4 interaction size
    unsigned int exclSz = app->topology->exclusions.getTable().size();

    nonbonded = new OpenMM::NonbondedForce(sz,0);//exclSz);
    system->addForce(nonbonded);

    //normal interactions
    for (unsigned int i = 0; i < sz; i++){
      int type1 = app->topology->atoms[i].type;
      Real sigma = app->topology->atomTypes[type1].sigma;
      //topo->atomTypes[i].sigma14 = par.nonbondeds[bi].sigma14;
      Real epsilon = app->topology->atomTypes[type1].epsilon;
      //topo->atomTypes[i].epsilon14 = par.nonbondeds[bi].epsilon14;
      Real charge = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
      Real mass = app->topology->atoms[i].scaledMass;
      //
      Real c6 = 4.0 * epsilon * pow(sigma, 6.0) * Constant::KCAL_KJ * 1e-6;
      Real c12 = 4.0 * epsilon * pow(sigma, 12.0) * Constant::KCAL_KJ * 1e-12;

#ifdef DEBUG
      mFile << c6 << " " << c12 << " " << charge << " " << mass << std::endl;  
#endif

      if (c12 <= 0){
			  nonbonded->setParticleParameters(i, charge, 1.0, 0.0);
      }else{
			  nonbonded->setParticleParameters(i, charge, pow(c12/c6, (1.0/6.0)), c6*c6/(4.0*c12));
      }

      system->setParticleMass(i, mass);

    }

#ifdef DEBUG
    mFile << std::endl;
#endif

    //1-4 interactions	
	  std::vector<NBForce> mForces;

    for (unsigned int i = 0; i < exclSz; i++){
      if ( (app->topology->exclusions.getTable())[i].excl == EXCLUSION_MODIFIED) {
        unsigned int atom1 = (app->topology->exclusions.getTable())[i].a1;
        unsigned int atom2 = (app->topology->exclusions.getTable())[i].a2;

        unsigned int type1 = app->topology->atoms[atom1].type;
        unsigned int type2 = app->topology->atoms[atom2].type;
        Real sigma = 0.5 * (app->topology->atomTypes[type1].sigma +
                              app->topology->atomTypes[type2].sigma);
        Real epsilon = sqrt(app->topology->atomTypes[type1].epsilon * 
                              app->topology->atomTypes[type2].epsilon);//0.5
        Real chargeij =  FudgeQQ * //app->topology->coulombScalingFactor *FudgeQQ
                          (app->topology->atoms[atom1].scaledCharge / Constant::SQRTCOULOMBCONSTANT) *
                            (app->topology->atoms[atom2].scaledCharge / Constant::SQRTCOULOMBCONSTANT); 

        Real c6 =  FudgeLJ * (4.0 * epsilon * pow(sigma, 6.0) * Constant::KCAL_KJ * 1e-6); //FudgeLJ
        Real c12 = FudgeLJ * (4.0 * epsilon * pow(sigma, 12.0) * Constant::KCAL_KJ * 1e-12);

        Real epsilon2 = (c6*c6)/(4.0*c12);
        Real sigma2 = pow(c12/c6,  (1.0/6.0));
		
		    mForces.push_back( NBForce( atom1, atom2, chargeij, sigma2, epsilon2, c6, c12 ) );

        //mFile << i << " " << atom1 << " " << atom2 << " " << chargeij << " " << 
        //    sigma2 << " " << epsilon2 << " " << c6 << " " << c12 << std::endl;  

        if (c12 <= 0) {
          //nonbonded->setNonbonded14Parameters(i, atom1, atom2, chargeij , 1.0, 0.0);
        } else {
          //nonbonded->setNonbonded14Parameters(i, atom1, atom2, chargeij, sigma2, epsilon2);
        }
      }
    }
	
	std::sort( mForces.begin(), mForces.end() );

#ifdef DEBUG
	mFile << "NonBonded 14 Force" << std::endl;

	for( unsigned int i = 0; i < mForces.size(); i++){
		const NBForce &temp = mForces[i];
		
		mFile  << temp.atom1 << " " << temp.atom2 << " " << temp.charge << " " << temp.sigma << " " << temp.epsilon << " " << temp.c6 << " " << temp.c12 << std::endl;
	}
	mFile << std::endl;
#endif

  }


  //openMM Initialize
  if( myIntegratorType == 1) {
    integrator = new OpenMM::LangevinIntegrator(myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS); //ps
  } else {
    //integrator = new OpenMM::NMLIntegrator(myLangevinTemperature, myGamma, getTimestep() * Constant::FS_PS, &app->eigenInfo); //ps
  }
  context = new OpenMM::OpenMMContext(*system, *integrator);

  OpenMM::Vec3 openMMvecp, openMMvecv;
  for (unsigned int i = 0; i < sz; ++i){
    for ( int j = 0; j < 3; j++){
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

  mFile.close();
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

  for (unsigned int i = 0; i < sz; ++i){
   for (int j = 0; j < 3; j++){
     app->positions[i].c[j] = openMMpositions[i][j] * Constant::NM_ANGSTROM; //nm to A
     app->velocities[i].c[j] = openMMvelocities[i][j];// * Constant::NM_ANGSTROM * Constant::TIMEFACTOR; //nm/ps to A/fs?
     (*myForces)[i].c[j] = openMMforces[i][j] * Constant::INV_NM_ANGSTROM * Constant::KJ_KCAL; //KJ/nm to Kcal/A
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
  parameters.push_back(Parameter( "RBDihedralForce", Value( RBDihedralForce, ConstraintValueType::NoConstraints() ), false ));
  parameters.push_back(Parameter( "PeriodicTorsion", Value( PeriodicTorsion, ConstraintValueType::NoConstraints() ), false ));
  parameters.push_back(Parameter( "NonbondedForce", Value( NonbondedForce, ConstraintValueType::NoConstraints() ), false ));
  parameters.push_back(Parameter( "IntegratorType", Value( myIntegratorType, ConstraintValueType::NotNegative() ), 1 ));

}

STSIntegrator *OpenMMIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  OpenMMIntegrator* myIntegP = new OpenMMIntegrator(values[0], fg);

  std::vector<Value> myValues(values.begin() + 1, values.end());

  myIntegP->setupValues(myValues);

  return (STSIntegrator*)myIntegP;

}

void OpenMMIntegrator::setupValues(std::vector<Value> &values) {

  //these must be in the same order as getParameters()
  myLangevinTemperature = values[0];
  myGamma = (Real)values[1] / (1000.0 * Constant::INV_TIMEFACTOR);
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed = values[2];
  HarmonicBondForce = values[3];
  HarmonicAngleForce = values[4];
  RBDihedralForce = values[5];
  PeriodicTorsion = values[6];
  NonbondedForce = values[7];
  myIntegratorType = values[8];

}

