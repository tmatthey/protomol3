#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/BlockMatrix.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include<iostream>
#include<fstream>
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
  //__________________________________________________ NormalModeDiagonalize

  const string NormalModeDiagonalize::keyword( "NormalModeDiagonalize" );

  NormalModeDiagonalize::NormalModeDiagonalize() : MTSIntegrator(), NormalModeUtilities(), eigAlloc(false), hessianCounter(0), rediagCounter(0)
  {
        eigAlloc=false;

  }

  NormalModeDiagonalize::NormalModeDiagonalize(int cycles, int redi,  bool fDiag, bool rRand, 
                                               Real redhy, Real eTh, int bvc, int rpb, Real dTh,
                                               ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), NormalModeUtilities( 1, 1, 91.0, 1234, 300.0), eigAlloc(false), fullDiag(fDiag), removeRand(rRand),
        rediagCount(redi), rediagHysteresis(redhy), 
        hessianCounter(0), rediagCounter(0), eigenValueThresh(eTh), blockVectorCols(bvc), residuesPerBlock(rpb), blockCutoffDistance(dTh)   //
  {
        eigAlloc=false;
        //
        rHsn.findForces(overloadedForces);	//find forces and parameters
  }

  NormalModeDiagonalize::~NormalModeDiagonalize() 
  {   
    //output stats
    report.precision(5);
    if(rediagCounter && hessianCounter){
        report <<plain<<"NML Timing: Hessian: "<<(blockDiag.hessianTime.getTime()).getRealTime()
            <<"[s] ("<<hessianCounter<<" times), diagonalize: "<<(blockDiag.rediagTime.getTime()).getRealTime()<<
                "[s] ("<<rediagCounter<<" re-diagonalizations)."<<endl;
        if(!fullDiag) report <<plain<<"NML Memory: Hessian: "<<memory_Hessian
            <<"[Mb], diagonalize: "<<memory_eigenvector<<
            "[Mb], vectors: "<<_3N*_rfM*sizeof(double)/1000000<<"[Mb]."<<endl;
    }
    //de-allocate
    if(mhQu!=NULL && eigAlloc) delete [] mhQu;
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
    //Initialize Hessian array, OR assign hessian array for residues.
    if(fullDiag){
      rHsn.initialData(_3N);
    }else{
      //assign hessian array for residues, and clear.
      bool fullE = false;
      if(blockCutoffDistance == 0.0) fullE = true;
      rHsn.initialResidueData(app->topology, residuesPerBlock, fullE); //added residuess per block
    }
    //should check to see if mhQ is large enough and re-use, then fix de-allocation
    if(mhQu != NULL){
        firstDiag = eigAlloc = false;
        validMaxEigv = true;
    }else{
        firstDiag = eigAlloc = true;
        validMaxEigv = false;
    }
    //assign arrays
    try{
        //should check to see if mhQ is large enough and re-use, then fix de-allocation
        if(mhQu == NULL){
            if(fullDiag) mhQu = new double[_3N*_3N];
            else mhQu = new double[_3N*_rfM];
        }
        //
    }catch(bad_alloc&){
        report << error << "Eigenvector array allocation error."<<endr;
    }
    //Initialize BlockHessianDiagonalize, pass BlockHessian if Blocks (not full diag)
    if(fullDiag){
      blockDiag.initialize(_3N);
    }else{
      blockDiag.initialize(&rHsn, _3N);
    }
    //Diagnostics
    memory_Hessian = memory_eigenvector = 0;
    //setup rediag counter in case valid
    nextRediag = 0;	//rediag first time 
    //save positions where diagonalized for checkpoint save (assume I.C. if file)
    diagAt = *&app->positions;
    //timers/counters for diagnostics
    hessianCounter = rediagCounter = rediagUpdateCounter = 0;
    //
  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeDiagonalize::run(int numTimesteps) {

    if( numTimesteps < 1 )
      return;
    //main loop
    app->energies.clear();
    for(int i=0;i<numTimesteps;i++){
        //Full diagonalization
        if((rediagCount && (int)(app->topology->time/getTimestep()) >= nextRediag) || firstDiag){
            nextRediag += rediagCount;
            //
            report <<debug(1)<<"[NormalModeDiagonalize::run] Finding diagonalized Hessian."<<endr;
            //save positions where diagonalized for checkpoint save
            diagAt = *&app->positions;
            //remove last random perturbation?
            if(removeRand) (*&app->positions).intoSubtract(myLastNormalMode->gaussRandCoord1);
            //Diagonalize
            if(fullDiag){
                //****Full method**********************************************************************//
                // Uses BLAS/LAPACK to do 'brute force' diagonalization                                //
                //*************************************************************************************//
                report << hint << "Start diagonalization."<<endr;
                //Find Hessians
                blockDiag.hessianTime.start();	//time Hessian
                rHsn.clear();
                rHsn.evaluate(&app->positions, app->topology, true);	//mass re-weighted hessian
                report << hint << "Hessian found."<<endr;
                //
                blockDiag.hessianTime.stop();	//stop timer
                hessianCounter++;
                //Diagonalize        
                blockDiag.rediagTime.start();
                int numeFound;                
                int info = blockDiag.diagHessian(mhQu, blockDiag.eigVal, rHsn.hessM, _3N, numeFound);
                if(info) report << error << "Full diagonalization failed."<<endr;
                //find number of -ve eigs
                int ii;
                for(ii=0;ii<_3N-3;ii++) if(blockDiag.eigVal[ii+3] > 0) break; 
                //
                report <<debug(1)<<"[NormalModeDiagonalize::run] Full diagonalize. No. negative eigenvales = "<<ii<<endr;
                for(int i=0;i<_3N;i++) blockDiag.eigIndx[i] = i;
                blockDiag.absSort(mhQu, blockDiag.eigVal, blockDiag.eigIndx, _3N);
                //
                blockDiag.rediagTime.stop();
                rediagCounter++;
                //set flags if firstDiag
                if(firstDiag){
                    numEigvectsu = _3N;
                    maxEigvalu = blockDiag.eigVal[_3N-1];
                    validMaxEigv = true;
                    firstDiag = false;
                }
            }else{
                //****Coarse method**************************************************************************//
                // Process:  Finds isolated 'minimized' block (of residues) Hessians   [evaluateResidues]    //
                //           Diagonalizes blocks to form block eigenvectors B          [findCoarseBlockEigs] //
                //           Finds actual Hessian H (but coarse grained) then S=B^THB  [innerHessian]        //
                //           Diagonalizes S to get eigenvectors Q, then approximate                          //
                //           eigenvectors are the first 'm' columns of BQ.                                   //
                //*******************************************************************************************//
                report << hint << "Start coarse diagonalization."<<endr;  
                Real max_eigenvalue = blockDiag.findEigenvectors(&app->positions, app->topology, 
                                                                  mhQu, _3N, _rfM, 
                                                                  blockCutoffDistance, eigenValueThresh, blockVectorCols);
                //Stats/diagnostics
                rediagCounter++; hessianCounter++;
                memory_Hessian = (rHsn.memory_base + rHsn.memory_blocks) * sizeof(Real) / 1000000;
                memory_eigenvector = blockDiag.memory_footprint * sizeof(Real) / 1000000;
                //set flags if firstDiag (firstDiag can now be coarse)
                if(firstDiag){
                    numEigvectsu = _3N;
                    //####FIX LEFT LIKE THIS FOR REGRESSION TESTS
                    maxEigvalu = 1000;//max_eigenvalue;	//maximum from blocks
                    validMaxEigv = true;
                    firstDiag = false;
                }
                //
                report << hint << "Coarse diagonalization complete. Maximum eigenvalue = " << max_eigenvalue << "." << endr;
            }
            //fix positions, Removed Rand? 
            if(removeRand)
                    *&app->positions = diagAt;	//back to original positions
            //sift current velocities/forces
            myNextNormalMode->subSpaceSift(&app->velocities, myForces);

        }
        //run integrator
        myNextIntegrator->run(1);//numTimesteps); //SHOULD DO 'rediagCount' STEPS HERE

    }

  } 

  //********************************************************************************************************************************************

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeDiagonalize::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("reDiagFrequency", Value(rediagCount,ConstraintValueType::NotNegative()),0,Text("Frequency of re-diagonalization (steps).")));
    parameters.push_back(Parameter("fullDiag",Value(fullDiag,ConstraintValueType::NoConstraints()),false,Text("Full diagonalization?")));
    parameters.push_back(Parameter("removeRand",Value(removeRand,ConstraintValueType::NoConstraints()),false,Text("Remove last random perturbation?")));
    parameters.push_back(Parameter("rediagHysteresis",Value(rediagHysteresis,ConstraintValueType::NotNegative()),0.0,Text("Re-diagonalization hysteresis.")));
    parameters.push_back(Parameter("eigenValueThresh",Value(eigenValueThresh,ConstraintValueType::NotNegative()),5.0,Text("'Inner' eigenvalue inclusion threshold.")));
    parameters.push_back(Parameter("blockVectorCols",Value(blockVectorCols,ConstraintValueType::NotNegative()),0,Text("Target number of block eigenvector columns.")));
    parameters.push_back(Parameter("residuesPerBlock", Value(residuesPerBlock,ConstraintValueType::NotNegative()),1,Text("Residues per block.")));
    parameters.push_back(Parameter("blockCutoffDistance", Value(blockCutoffDistance,ConstraintValueType::NotNegative()),10,Text("Block cutoff distance for electrostatic forces.")));
  }

  MTSIntegrator* NormalModeDiagonalize::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeDiagonalize(values[0],values[1],values[2],values[3],values[4],values[5],values[6],
      values[7], values[8], fg,nextIntegrator);
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
