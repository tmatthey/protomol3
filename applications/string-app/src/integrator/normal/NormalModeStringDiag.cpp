#include <src/integrator/normal/NormalModeStringDiag.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <iostream>
#include <fstream>
#include <iomanip>
//#include "ModifierForceProjection.h"

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#else
#if defined (HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeStringDiag

  const string NormalModeStringDiag::keyword( "NormalModeStringDiag" );

  NormalModeStringDiag::NormalModeStringDiag() : MTSIntegrator(), NormalModeUtilities(), eigAlloc(false), hessianCounter(0), rediagCounter(0)
  {
        eigVal=NULL;eigIndx=NULL;innerEigVec=NULL;innerEigVal=NULL;innerHess=NULL;eigAlloc=false;
        T1=NULL;HQ=NULL;tempMxM=NULL;temp3NxM=NULL;w=NULL;
        orig_Eig_Vec = NULL;
  }

  NormalModeStringDiag::NormalModeStringDiag(int cycles, int avs,  Real avss, int redi, bool fDiag, bool rRand, 
                                                    int mins, Real minl, Real redn, Real redhy, Real spd, int maxi, bool rBond,int doM,
                                                        ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), NormalModeUtilities( 1, 1, 91.0, 1234, 300.0), eigAlloc(false), fullDiag(fDiag), removeRand(rRand), noAvStep(avs), 
        avStep(avss), rediagCount(redi), minSteps(mins), minLim(minl), 
      rediagThresh(redn), rediagHyst(redhy), spdOff(spd), maxIterations(maxi), hessianCounter(0), rediagCounter(0), removeBondedEigs(rBond), doMin(doM)   //
  {
        eigVal=NULL;eigIndx=NULL;innerEigVec=NULL;innerEigVal=NULL;innerHess=NULL;eigAlloc=false;
        T1=NULL;HQ=NULL;tempMxM=NULL;temp3NxM=NULL;w=NULL;
        orig_Eig_Vec = NULL;
        //
        hsn.findForces(overloadedForces);	//find forces and parameters
  }


  NormalModeStringDiag::~NormalModeStringDiag() 
  {   
    //output stats
    report.precision(5);
    if(rediagCounter && hessianCounter){
        rediagIters /= rediagCounter;
        report <<plain<<"NML Timing: Hessian: "<<(hessianTime.getTime()).getRealTime()
            <<"[s] ("<<hessianCounter<<" times), diagonalize: "<<(rediagTime.getTime()).getRealTime()<<
                "[s] ("<<rediagCounter<<" tests, "<<rediagUpdateCounter<<" re-diag, "<<rediagIters<<" iterations)."<<endl;
    }
    //de-allocate
    if(mhQu!=NULL && eigAlloc) delete [] mhQu;
    if(eigVal!=NULL) delete [] eigVal;
    if(eigIndx!=NULL) delete [] eigIndx;
    if(innerEigVec!=NULL) delete [] innerEigVec;
    if(innerEigVal!=NULL) delete [] innerEigVal;
    if(innerHess!=NULL) delete [] innerHess;
    if(T1!=NULL) delete [] T1;
    if(HQ!=NULL) delete [] HQ;
    if(tempMxM!=NULL) delete [] tempMxM;
    if(temp3NxM!=NULL) delete [] temp3NxM;
    if(w!=NULL) delete [] w;

    if (orig_Eig_Vec != NULL) delete [] orig_Eig_Vec;
  }

  void NormalModeStringDiag::initialize(ProtoMolApp *app) {
    MTSIntegrator::initialize(app);
    //
    myNextNormalMode  = dynamic_cast<NormalModeUtilities*>(myNextIntegrator);
    //Using complement of next integrator, so copy 
    firstMode = myNextNormalMode->firstMode; numMode = myNextNormalMode->numMode;
    //NM initialization, but fix dof, as set by next integrator
    NormalModeUtilities::initialize((int)app->positions.size(), app->topology, myForces, COMPLIMENT_FORCES); //last for complimentary forces, no gen noise
    _rfM = myNextNormalMode->_rfM;
    _idM = _rfM + 2; /* only used for traceMin */
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
    innerEigVec = new double[_idM*_idM];
    innerEigVal = new double[_idM];
    innerHess = new double[_idM*_idM];
    //steepest decent calculation variables
    T1 = new double[_3N*_idM];
    HQ = new double[_3N*_idM];
    tempMxM = new double[_idM*_idM];
    temp3NxM = new double[_3N*_idM];
    w = new double[_idM];
    //
    if(mhQu==NULL || eigVal==NULL || eigIndx==NULL || hsn.hessM==NULL || innerEigVec==NULL || 
        innerEigVal==NULL || innerHess==NULL || T1==NULL || HQ==NULL || tempMxM==NULL || temp3NxM==NULL || w==NULL) 
            report << error << "Eigenvector array allocation error."<<endr;

    orig_Eig_Vec = new double[_3N*_3N];
    for(int i=0;i<_rfM*_3N;i++) orig_Eig_Vec[i] = mhQu[i];
    //
    //setup rediag counter incase valid
    nextRediag = 0;	//rediag first time 
    //
    if(rediagHyst > rediagThresh) rediagHyst = rediagThresh;
    //save positions where diagonalized for checkpoint save (assume I.C. if file)
    diagAt = *&app->positions;
    myAveragePos = *&app->positions;
    numAveragePos = 1;
    //timers/counters for diagnostics
    rediagTime.reset();
    hessianTime.reset();
    hessianCounter = rediagCounter = rediagIters = rediagUpdateCounter = 0;


  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeStringDiag::run(int numTimesteps) {
    if( numTimesteps < 1 )
      return;
    //main loop
    app->energies.clear();
    for(int i=0;i<numTimesteps;i++){
        myNextIntegrator->run(myCycleLength);

        if (doMin) minimizer(minLim, minSteps, true, false, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);

    }

  }  


//********************************************************************************************************************************************


  void NormalModeStringDiag::removeBondEig(){
    int a1, a2, offs1, offs2;
    Vector3D aDiff;
    Real dotProd;

    for (unsigned int i=0; i<app->topology->bonds.size(); i++){
        //find position vectors and find the normed difference
        a1 = app->topology->bonds[i].atom1; a2 = app->topology->bonds[i].atom2;
        if((app->topology->atoms[a1].name.c_str())[0] == 'H' || (app->topology->atoms[a2].name.c_str())[0] == 'H'){
            //report << hint << "[NormalModeStringDiag] i= "<<i<<" a1= "<<a1<<" a2 = "<<a2<<" name1= "<<app->topology->atoms[a1].name<<" name2= "<<app->topology->atoms[a2].name<<endr;
            aDiff = app->positions[a2]-app->positions[a1];
            //find ratio of masses
            Real rMass = app->topology->atoms[a1].scaledMass / app->topology->atoms[a2].scaledMass;
            Real norm1 = sqrt(1.0/(1.0+rMass));
            Real norm2 = sqrt(1.0/(1.0+1.0/rMass));
            Real aDnorm = aDiff.norm();
            if(aDnorm) aDiff /= aDnorm;
            //remove from eigenvector set 
            for(int j=0; j<_idM; j++){
                //point to first double in array for both atoms
                offs1 = j*_3N + a1*3;
                offs2 = j*_3N + a2*3;
                //find dot product
                dotProd = 0.0;
                for(int k=0;k<3;k++){
                    dotProd += mhQu[offs1+k]*aDiff.c[k]*norm1;
                    dotProd -= mhQu[offs2+k]*aDiff.c[k]*norm2;
                }
                //remove dotprodct times eigvectors
                for(int k=0;k<3;k++){
                    mhQu[offs1+k] -= dotProd*aDiff.c[k]*norm1;
                    mhQu[offs2+k] += dotProd*aDiff.c[k]*norm2;
                }
                //report << hint << "[NormalModeStringDiag] dp1 = "<<dotProd<<" dp2 = "<<dotProd2<<endr;
            }
        }
        //Real r_0 = myTopo->bonds[i].restLength;
        //Real k = myTopo->bonds[i].springConstant;
    }
    //Normalize the resulting Q
    for(int j=0;j<_idM;j++){
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
        int incxy = 1;
        int n = _3N;
#endif
        double norm = 0;
#if defined(HAVE_LAPACK)   
        norm = dnrm2_(&n, &mhQu[j*_3N], &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
        norm = dnrm2_(n, &mhQu[j*_3N], incxy);
#endif
#endif
        if(norm == 0.0) report << error << "[NormalModeStringDiag] Vector with zero norm."<<endr;
        for(int i=0;i<_3N;i++) mhQu[i + j*_3N] /= norm; 
    }
  }


//********************************************************************************************************************************************

  //NVE setps for Hessian average
  void NormalModeStringDiag::VerletAverage(){
    Real h = avStep * Constant::INV_TIMEFACTOR;

    myAveragePos = *&app->positions;
    numAveragePos = 1;
    calculateForces();
    //Verlet prop.
    for(int nos=1;nos<noAvStep;nos++){
        report <<debug(5)<<"[NormalModeStringDiag::run] averaging step = "<<nos<<" step size = "<<avStep<<endr;
        for (int i = 0; i < _N; ++i ) {
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
          app->positions[i] += app->velocities[i] * h;
        } 
        for(int j=0;j<_N;j++) myAveragePos[j] += app->positions[j];
        numAveragePos++;
        //hsn.evaluate(myPositions, myTopo, true);	//true for mass re-weight;
        calculateForces();
        for (int i = 0; i < _N; ++i ) {
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
        } 
    }
    //Average 
    for(int j=0;j<_N;j++) myAveragePos[j] /= (Real)numAveragePos;
    //for(int i=0 ; i<_3N*_3N ; i++) hsn.hessM[i] /= (double)noAvStep;	//divide sum of Hessians

  }

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeStringDiag::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("averageSteps", Value(noAvStep,ConstraintValueType::NotNegative()),1,Text("Hessian averaged over number of steps.")));
    parameters.push_back(Parameter("avStepSize",Value(avStep,ConstraintValueType::NotNegative()),1.0,Text("Step size for Hessian averaging.")));
    parameters.push_back(Parameter("reDiagTestFrequency", Value(rediagCount,ConstraintValueType::NotNegative()),0,Text("Frequency of re-diagonalization (steps).")));
    parameters.push_back(Parameter("fullDiag",Value(fullDiag,ConstraintValueType::NoConstraints()),false,Text("Full diagonalization?")));
    parameters.push_back(Parameter("removeRand",Value(removeRand,ConstraintValueType::NoConstraints()),false,Text("Remove last random perturbation?")));
    parameters.push_back(Parameter("minSteps", Value(minSteps,ConstraintValueType::NotNegative()),0,Text("Max. number of minimizer steps.")));
    parameters.push_back(Parameter("minLim",Value(minLim,ConstraintValueType::NotNegative()),1.0,Text("Minimization limit kcal mol^{-1}.")));
    parameters.push_back(Parameter("upperThreshold",Value(rediagThresh,ConstraintValueType::NotNegative()),0.0,Text("Re-diagonaliztion threshold.")));
    parameters.push_back(Parameter("lowerThreshold",Value(rediagHyst,ConstraintValueType::NotNegative()),0.0,Text("Re-diagonaliztion threshold with hysteresis.")));
    parameters.push_back(Parameter("spdOffset",Value(spdOff,ConstraintValueType::NotNegative()),200.0,Text("Factor to ensure Hessisn SPD.")));
    parameters.push_back(Parameter("maxIterations", Value(maxIterations,ConstraintValueType::NotNegative()),1000,Text("Maximum re-diagonalization iterations.")));
    parameters.push_back(Parameter("removeBondedEigs",Value(removeBondedEigs,ConstraintValueType::NoConstraints()),false,Text("Remove bonded eigenvectors?")));
    parameters.push_back(Parameter("domin",Value(doMin,ConstraintValueType::NotNegative()),0,Text("Minimize before reparameterize?")));
  }

  MTSIntegrator* NormalModeStringDiag::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeStringDiag(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],
        values[8],values[9],values[10],values[11],values[12],values[13],fg,nextIntegrator);
  }

  //void NormalModeStringDiag::addModifierAfterInitialize(){
  //  adoptPostForceModifier(new ModifierForceProjection(this));
  //  MTSIntegrator::addModifierAfterInitialize();
  //}

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************
  void NormalModeStringDiag::utilityCalculateForces(){
      app->energies.clear();

      calculateForces();
  }

}
