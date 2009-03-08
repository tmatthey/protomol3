#include <src/integrator/normal/NormalModeFERelax.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>


using namespace std;


using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeFERelax

  const string NormalModeFERelax::keyword( "NormalModeFERelax" );

  NormalModeFERelax::NormalModeFERelax() : MTSIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeFERelax::NormalModeFERelax(int cycles, Real minimlim, bool rediag, bool simplemin, bool dm,
                                    ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator)
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
        NormalModeUtilities( 1, 1, 91.0, 1234, 300.0),
            minLim(minimlim), reDiag(rediag), simpleMin(simplemin), doMin(dm)
  {
  }


  NormalModeFERelax::~NormalModeFERelax() 
  {  
  }

  void NormalModeFERelax::initialize(ProtoMolApp* app){
    MTSIntegrator::initialize(app);
    //test not topmost integrator
    if(top() == this) report << error << "NormalModeFERelax cannot be top integrator."<<endr;
    //
    initializeForces();
    //
    myPreviousNormalMode  = dynamic_cast<NormalModeUtilities*>(myPreviousIntegrator);
    //check valid eigenvectors
    firstMode = myPreviousNormalMode->firstMode; numMode = myPreviousNormalMode->numMode;
    //NM initialization if OK
    NormalModeUtilities::initialize((int)app->positions.size(), app->topology, myForces, COMPLIMENT_FORCES); //last for complimentary forces
    //
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    //
    //diagnostics
    avItrs = 0;		//average number of minimizer iterations/force calcs
    avMinForceCalc = 0;
    numSteps = 0;	//total steps

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeFERelax::run(int numTimesteps) {

    if( numTimesteps < 1 )
      return;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    //do minimization with local forces, max loop 100, set subSpace minimization true
   
    Real minPotEnergy; 
    if (doMin) {
    //cout<<"Minimization in relax..."<<endl;
    itrs = minimizer(minLim, 100, simpleMin, reDiag, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);
     minPotEnergy = app->energies.potentialEnergy();//save potential energy before random
   }
    //constraints?
    app->energies.clear();
    //run brownian if any remaining modes
    save_pos = *&app->positions;
    if(testRemainingModes() && (myCycleLength>1)) myNextIntegrator->run(myCycleLength);	//cyclelength 
    else cout<<"Not doing brownian "<<endl;
    *&app->positions =  save_pos;
    if (doMin) (app->energies)[ScalarStructure::OTHER] = minPotEnergy;	//store minimum potential energy
    if(*Q == NULL){	//rediagonalize?
        if(myPreviousIntegrator == NULL) 
            report << error << "[NormalModeFERelax::Run] Re-diagonalization forced with NormalModeFERelax as outermost Integrator. Aborting."<<endr;
        return;
    }
    //
    postStepModify();
    //
  }  

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeFERelax::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));
    parameters.push_back(Parameter("rediag",Value(reDiag,ConstraintValueType::NoConstraints()),false,Text("Force re-digonalize")));
    parameters.push_back(Parameter("simplemin",Value(simpleMin,ConstraintValueType::NoConstraints()),true,Text("Simple minimizer or exact minima projection.")));
    parameters.push_back(Parameter("doMinimize",Value(doMin,ConstraintValueType::NoConstraints()),true,Text("Perform minimization before Brownian?")));
  }

  MTSIntegrator* NormalModeFERelax::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeFERelax(values[0],values[1],values[2],values[3],values[4],fg,nextIntegrator);
  }

  //void NormalModeFERelax::addModifierAfterInitialize(){
  //  adoptPostForceModifier(new ModifierForceProjection(this));
  //  MTSIntegrator::addModifierAfterInitialize();
  //}

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************

  void NormalModeFERelax::utilityCalculateForces(){
      //cout<<"Inside utilityCalculateForces of NormalModerelax "<<endl;
      app->energies.clear();	//need this as MTS, not innermost
      calculateForces();
  }

}
