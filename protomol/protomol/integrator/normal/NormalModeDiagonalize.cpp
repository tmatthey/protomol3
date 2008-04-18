#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
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
  //__________________________________________________ NormalModeDiagonalize

  const string NormalModeDiagonalize::keyword( "NormalModeDiagonalize" );

  NormalModeDiagonalize::NormalModeDiagonalize() : MTSIntegrator(), NormalModeUtilities()
  {
      eigVal=NULL;eigIndx=NULL;eigAlloc=false;origEigVec=NULL;
  }

  NormalModeDiagonalize::NormalModeDiagonalize(int cycles, int avs,  Real avss, int redi, bool fDiag, bool rRand,
                                                    int mins, Real minl, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), NormalModeUtilities( 1, 1, 91.0, 1234, 300.0), fullDiag(fDiag),  
        removeRand(rRand), noAvStep(avs), avStep(avss), rediagCount(redi), minSteps(mins), minLim(minl)    //
  {
        eigVal=NULL;eigIndx=NULL;innerEigVec=NULL;innerEigVal=NULL;innerHess=NULL;eigAlloc=false;origEigVec=NULL;
        //
        hsn.findForces(overloadedForces);	//find forces and parameters
  }


  NormalModeDiagonalize::~NormalModeDiagonalize() 
  {  
    if(mhQu!=NULL && eigAlloc) delete [] mhQu;
    if(eigVal!=NULL) delete [] eigVal;
    if(eigIndx!=NULL) delete [] eigIndx;
    if(innerEigVec!=NULL) delete [] innerEigVec;
    if(innerEigVal!=NULL) delete [] innerEigVal;
    if(innerHess!=NULL) delete [] innerHess;
    if(origEigVec!=NULL) delete [] origEigVec;
  }

  void NormalModeDiagonalize::initialize(ProtoMolApp *app) {
    MTSIntegrator::initialize(app);
    //test topmost integrator
    if(top() != this) report << error << "NormalModeDiagonalize not top integrator."<<endr;
    //
    myNextNormalMode  = dynamic_cast<NormalModeUtilities*>(myNextIntegrator);
    //Using complement of next integrator, so copy 
    firstMode = myNextNormalMode->firstMode; numMode = myNextNormalMode->numMode;
    //NM initialization, but fix dof, as set by next integrator
    NormalModeUtilities::initialize((int)app->positions.size(), app->topology, myForces, COMPLIMENT_FORCES); //last for complimentary forces, no gen noise
    _rfM = myNextNormalMode->_rfM;
    //
    myLastNormalMode  = dynamic_cast<NormalModeUtilities*>(bottom());
    //do first force calculation, and remove non sub-space part
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    initializeForces();
    //
    //should check to see if mhQ is large enough and re-use, then fix de-allocation
    if(mhQu != NULL){
        firstDiag = eigAlloc = false;
        validMaxEigv = true;
    }else{
        mhQu = new double[_3N*_3N];
        firstDiag = eigAlloc = true;
        validMaxEigv = false;
    }
    //
    eigVal = new double[_3N];
    eigIndx = new int[_3N];
    //initialize Hessian array
    hsn.initialData(_3N);
    //initialise inner Hessian arrays
    innerEigVec = new double[_rfM*_rfM];
    innerEigVal = new double[_rfM];
    innerHess = new double[_rfM*_rfM];
    //
    origEigVec = new double[_rfM*_3N];
    //
    if(mhQu==NULL || eigVal==NULL || eigIndx==NULL || hsn.hessM==NULL) report << error << "Eigenvector array allocation error."<<endr;
    //
    for(int i=0;i<_rfM*_3N;i++) origEigVec[i] = mhQu[i];
    //setup rediag counter incase valid
    nextRediag = 0;//rediagCount;	//rediad first time 
    //save positions where diagonalized for checkpoint save (assume I.C. if file)
    diagAt = *&app->positions;

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeDiagonalize::run(int numTimesteps) {
    if( numTimesteps < 1 )
      return;
    if((rediagCount && (int)(app->topology->time/getTimestep()) >= nextRediag) || firstDiag){
        if(!firstDiag) nextRediag += rediagCount;
        //Diagonalize if no input file
        report <<debug(1)<<"[NormalModeDiagonalize::run] Finding diagonalized Hessian."<<endr;
        //save positions where diagonalized for checkpoint save
        diagAt = *&app->positions;
        //index for sort by absolute
        for(int i=0;i<_3N;i++) eigIndx[i] = i;
        //mass re-weighted hessian
        hsn.clear();
        //NOTE: cannot minimize on first diagonalization as no max eigenvalue.
        //      Structure MUST be minimized.
        //remove last random perturbation?
        if(removeRand) app->positions.intoSubtract(myLastNormalMode->gaussRandCoord1);
        //minimizer AND valid maximum eigenvalue
	if(minSteps > 0 && validMaxEigv) {
	minimizer(minLim, minSteps, true, false, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);}
        //
        hsn.evaluate(&app->positions, app->topology, true);	//true for mass re-weight;
        //Avergaged?
        if(noAvStep > 1){
            Vector3DBlock tmpVel = *&app->velocities;	//save velocities
            VerletAverage();
            *&app->velocities = tmpVel;
        }
        //minimized/averaged/removed Rand? then reset positions
        if((minSteps > 0 && validMaxEigv) || noAvStep > 1 || removeRand)
                *&app->positions = diagAt;	//back to original positions
        //diagonalize
        if(!firstDiag && !fullDiag){	//inner diag
            getInnerHess(origEigVec, hsn.hessM, innerHess);	//find inner Hessian
            int info = diagHessian(innerEigVec, innerEigVal, innerHess, _rfM);	//diagonalize inner Hessian
            if( info == 0 ){
                //find max eigenvalue
                double maxe = 0.0;
                for(int l=0;l<_rfM;l++) if(innerEigVal[l]>maxe) maxe = innerEigVal[l];
                //find number of -ve eigs
                int ii;
                for(ii=0;ii<_3N-3;ii++) if(innerEigVal[ii+3] > 0) break; 
                //
                report <<debug(1)<<"[NormalModeDiagonalize::run] Inner diagonalized. Neg. eigenvalues = "<<ii<<" Maximum eigenvalue = "<<maxe<<endr;
                //absSort(innerEigVec, innerEigVal, eigIndx, _rfM);	//Don't need this
            }else{
                report << error << "Inner diagonalization failed."<<endr;
            }
            getNewEigs(mhQu, origEigVec, innerEigVec);	//calculate new eigenVectors
        }else{	//full diag
            int info = diagHessian(mhQu, eigVal, hsn.hessM, _3N);
            if( info == 0 ){
                //find number of -ve eigs
                int ii;
                for(ii=0;ii<_3N-3;ii++) if(eigVal[ii+3] > 0) break; 
                //
                report <<debug(1)<<"[NormalModeDiagonalize::run] Full diagonalized. No. negative eigenvales = "<<ii<<endr;
                absSort(mhQu, eigVal, eigIndx, _3N);
            }else{
                report << error << "Full diagonalization failed."<<endr;
            }
            numEigvectsu = _3N;
            maxEigvalu = eigVal[_3N-1];
            validMaxEigv = true;
            firstDiag = false;
        }
        //sift current velocities/forces
        myNextNormalMode->subSpaceSift(&app->velocities, myForces);
    }
    //main loop
    app->energies.clear();
    myNextIntegrator->run(numTimesteps);

  }  

  //NVE setps for Hessian average
  void NormalModeDiagonalize::VerletAverage(){
    Real h = avStep * Constant::INV_TIMEFACTOR;

    calculateForces();
    //Verlet prop.
    for(int nos=1;nos<noAvStep;nos++){
        report <<debug(5)<<"[NormalModeDiagonalize::run] averaging step = "<<nos<<" step size = "<<avStep<<endr;
        for( unsigned int i = 0; i < (unsigned int )_N; ++i ) {
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
          app->positions[i] += app->velocities[i] * h;
        } 
        hsn.evaluate(&app->positions, app->topology, true);	//true for mass re-weight;
        calculateForces();
        for( unsigned int i = 0; i < (unsigned int )_N; ++i ) {
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
        } 
    }
    //Avwerage hess
    for(int i=0 ; i<_3N*_3N ; i++) hsn.hessM[i] /= (double)noAvStep;	//divide sum of Hessians

  }

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeDiagonalize::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("averageSteps", Value(noAvStep,ConstraintValueType::NotNegative()),1,Text("Hessian averaged over number of steps.")));
    parameters.push_back(Parameter("avStepSize",Value(avStep,ConstraintValueType::NotNegative()),1.0,Text("Step size for Hessian averaging.")));
    parameters.push_back(Parameter("reDiagFrequency", Value(rediagCount,ConstraintValueType::NotNegative()),0,Text("Frequency of re-diagonalization (steps).")));
    parameters.push_back(Parameter("fullDiag",Value(fullDiag,ConstraintValueType::NoConstraints()),false,Text("Full diagonalization?")));
    parameters.push_back(Parameter("removeRand",Value(removeRand,ConstraintValueType::NoConstraints()),false,Text("Remove last random perturbation?")));
    parameters.push_back(Parameter("minSteps", Value(minSteps,ConstraintValueType::NotNegative()),0,Text("Max. number of minimizer steps.")));
    parameters.push_back(Parameter("minLim",Value(minLim,ConstraintValueType::NotNegative()),1.0,Text("Minimization limit kcal mol^{-1}.")));
  }

  MTSIntegrator* NormalModeDiagonalize::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeDiagonalize(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],fg,nextIntegrator);
  }

  //void NormalModeDiagonalize::addModifierAfterInitialize(){
  //  adoptPostForceModifier(new ModifierForceProjection(this));
  //  MTSIntegrator::addModifierAfterInitialize();
  //}

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************
  void NormalModeDiagonalize::utilityCalculateForces(){
	  app->energies.clear();

      calculateForces();
  }

}
