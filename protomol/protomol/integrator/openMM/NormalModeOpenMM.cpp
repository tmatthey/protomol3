#include <protomol/integrator/openMM/NormalModeOpenMM.h>
#include <protomol/base/Report.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

//#include "ModifierForceProjection.h"

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;


namespace ProtoMol {
  //__________________________________________________ NormalModeOpenMM

  const string NormalModeOpenMM::keyword( "NormalModeOpenMM" );

  NormalModeOpenMM::NormalModeOpenMM() : OpenMMIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeOpenMM::NormalModeOpenMM(Real timestep, int firstmode, int nummode, Real gamma, Real temperature, 
      Real minimlim, ForceGroup *overloadedForces) 
    : OpenMMIntegrator(timestep,overloadedForces), NormalModeUtilities( firstmode, nummode, gamma, 1234, temperature),
      minLim(minimlim)
  {
  }

  NormalModeOpenMM::~NormalModeOpenMM() 
  {  

  }

  void NormalModeOpenMM::initialize(ProtoMolApp* app){
    OpenMMIntegrator::initialize(app);
    initializeForces();
    //NM initialization
    NormalModeUtilities::initialize((int)app->positions.size(), app->topology,
				    myForces, NO_NM_FLAGS);
    //

    //Set up minimum limit
    app->eigenInfo.myMinimumLimit = minLim;
    //

  }

  void NormalModeOpenMM::run(int numTimesteps) {
    if( numTimesteps < 1 )
        return;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    preStepModify();
    //
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

    //
    postStepModify();
  }  

  void NormalModeOpenMM::getParameters(vector<Parameter>& parameters) const {
    OpenMMIntegrator::getParameters(parameters);

    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(NormalModeUtilities::myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));

  }

  STSIntegrator* NormalModeOpenMM::doMake(const vector<Value>& values,ForceGroup* fg)const{

    unsigned int numPar = OpenMMIntegrator::getParameterSize();

    std::vector<Value> ommValues(values.begin() + 1, values.begin() + numPar);

    NormalModeOpenMM* myIntegP = new NormalModeOpenMM(values[0],values[7],values[8],values[9],values[10],values[11],fg);

    myIntegP->OpenMMIntegrator::setupValues(ommValues);
  
    return (STSIntegrator*)myIntegP;

  }

  unsigned int NormalModeOpenMM::getParameterSize(){

    return OpenMMIntegrator::getParameterSize() + 5;

  }


}

