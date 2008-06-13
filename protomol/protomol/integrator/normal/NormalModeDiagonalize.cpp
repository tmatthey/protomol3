#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
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

  NormalModeDiagonalize::NormalModeDiagonalize() : MTSIntegrator(), NormalModeUtilities()
  {
        eigVal=NULL;eigIndx=NULL;innerEigVec=NULL;innerEigVal=NULL;innerHess=NULL;eigAlloc=false;
        T1=NULL;HQ=NULL;tempMxM=NULL;temp3NxM=NULL;w=NULL;
  }

  NormalModeDiagonalize::NormalModeDiagonalize(int cycles, int avs,  Real avss, int redi, bool fDiag, bool rRand, 
                                                    int mins, Real minl, Real redn, Real redhy, Real spd, int maxi, bool rBond,
                                                        ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), NormalModeUtilities( 1, 1, 91.0, 1234, 300.0), fullDiag(fDiag), removeRand(rRand), noAvStep(avs), 
        avStep(avss), rediagCount(redi), minSteps(mins), minLim(minl), 
            rediagThresh(redn), rediagHyst(redhy), spdOff(spd), maxIterations(maxi), removeBondedEigs(rBond)   //
  {
        eigVal=NULL;eigIndx=NULL;innerEigVec=NULL;innerEigVal=NULL;innerHess=NULL;eigAlloc=false;
        T1=NULL;HQ=NULL;tempMxM=NULL;temp3NxM=NULL;w=NULL;
        //
        hsn.findForces(overloadedForces);	//find forces and parameters
  }


  NormalModeDiagonalize::~NormalModeDiagonalize() 
  {   
    //output stats
    report.precision(5);
    if(rediagCounter) rediagIters /= rediagCounter;
    report <<plain<<"NML Timing: Hessian: "<<(hessianTime.getTime()).getRealTime()
        <<"[s] ("<<hessianCounter<<"), diagonalize: "<<(rediagTime.getTime()).getRealTime()<<"[s] ("<<rediagCounter<<", "<<rediagIters<<")."<<endl;
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
    _idM = _rfM + 2;
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
    hessianCounter = rediagCounter = rediagIters = 0;
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
            //Diagonalize if no input file
            report <<debug(1)<<"[NormalModeDiagonalize::run] Finding diagonalized Hessian."<<endr;
            //save positions where diagonalized for checkpoint save
            diagAt = *&app->positions;
            //mass re-weighted hessian
            hsn.clear();
            //NOTE: cannot minimize on first diagonalization as no max eigenvalue.
            //      Structure MUST be minimized.
            //remove last random perturbation?
            if(removeRand) (*&app->positions).intoSubtract(myLastNormalMode->gaussRandCoord1);
            //minimizer AND valid maximum eigenvalue
            if(minSteps > 0 && validMaxEigv)
                minimizer(minLim, minSteps, true, false, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);
            //
            //Avergaged?
            //time Hessian and average
            hessianTime.start();
            if(noAvStep > 1){
                myAveragePos = *&app->positions;
                numAveragePos = 1;
                Vector3DBlock tmpVel = *&app->velocities;	//save velocities
                VerletAverage();
                *&app->velocities = tmpVel;
                hsn.evaluate(&myAveragePos, app->topology, true);	//true for mass re-weight;
            }else{
                hsn.evaluate(&app->positions, app->topology, true);
            }
            hessianTime.stop();	//stop timer
            hessianCounter++;
            //minimized/averaged/removed Rand? then reset positions
            if((minSteps > 0 && validMaxEigv) || noAvStep > 1 || removeRand)
                    *&app->positions = diagAt;	//back to original positions
            //diagonalize
            rediagTime.start();
            //do full if set or no eig file
            if(fullDiag || firstDiag){
                int numeFound;
                int info = diagHessian(mhQu, eigVal, hsn.hessM, _3N, numeFound);
                if(info) report << error << "Full diagonalization failed."<<endr;
                //find number of -ve eigs
                int ii;
                for(ii=0;ii<_3N-3;ii++) if(eigVal[ii+3] > 0) break; 
                //
                report <<debug(1)<<"[NormalModeDiagonalize::run] Full diagonalize. No. negative eigenvales = "<<ii<<endr;
                for(int i=0;i<_3N;i++) eigIndx[i] = i;
                absSort(mhQu, eigVal, eigIndx, _3N);
                //set flags if firstDiag
                if(firstDiag){
                    numEigvectsu = _3N;
                    maxEigvalu = eigVal[_3N-1];
                    validMaxEigv = true;
                    firstDiag = false;
                }
            }else{	//else do trace method
                int iters = traceReDiagonalize(maxIterations, spdOff, _idM, removeBondedEigs); //true=remove bond eigs
                rediagIters += iters;
                if(iters == 0){	//exceeded maximum iterations, then do full diag
                    int numeFound;
                    int info = diagHessian(mhQu, eigVal, hsn.hessM, _3N, numeFound);
                    if(info) report << error << "Full diagonalization failed after maximum itereations exceeded."<<endr;
                    else report << hint << "Full diagonalization completed after maximum itereations exceeded."<<endr;
                    for(int i=0;i<_3N;i++) eigIndx[i] = i;
                    absSort(mhQu, eigVal, eigIndx, _3N);
                }
            }
            rediagCounter++;
            rediagTime.stop();
            //
            //sift current velocities/forces
            myNextNormalMode->subSpaceSift(&app->velocities, myForces);
        }
        //run integrator
        myNextIntegrator->run(1);//numTimesteps); //SHOULD DO 'rediagCount' STEPS HERE
    }

  }  

//********************************************************************************************************************************************

  int NormalModeDiagonalize::traceReDiagonalize(int maxIter, double spd, int idM, bool remBondEig){
    double traceThreshold, lambdaMax = 0, lambdaMax2, lambdaMin, spdOffs;
    int iter;
    //Blas variables
    char *transA, *transB;
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
    const char *transC, *transD;
    int nrhs;
    int info;
#endif
    int m, n, k, lda, ldb, ldc;
    double alpha, beta, norm = 0;

    //Set threshold for e-diag
    traceThreshold = rediagThresh;
    //iterations
    iter = 0;
    //############ Make spd #################################################
    spdOffs = spd;
    for(int i=0;i<_3N;i++) hsn.hessM[i+i*_3N] += spdOffs;
    //############ Remove bond eigenvectors from Q ##########################
    //removeBondEig();
    //############ Loop until converges or exceeds max iterations ###########
    while(iter < maxIter){
        iter++;
        //######## Remove bond eigenvectors from Q ##########################
        if(remBondEig) removeBondEig(); 
        //######## Given Q_i in mhQu, find \bar{Q}_i \in \bar{Q}^* ##########
        //Find S_i \in R^{mxm} (innerHess)
        //find inner Hessian and save intermediate HQ
        transA = "N"; transB = "N"; m = _3N; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = _3N; alpha = 1.0; beta = 0.0;
#if defined(HAVE_LAPACK) 
        dgemm_ (transA, transB, &m, &n, &k, &alpha, hsn.hessM, &lda, mhQu, &ldb, &beta, HQ, &ldc);
        transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
        dgemm_ (transA, transB, &m, &n, &k, &alpha, mhQu, &lda, HQ, &ldb, &beta, innerHess, &ldc);
#else
#if defined(HAVE_SIMTK_LAPACK)
        int len_trans = 1;
        dgemm_ (*transA, *transB, m, n, k, alpha, hsn.hessM, lda, mhQu, ldb, beta, HQ, ldc, len_trans, len_trans);
        transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
        dgemm_ (*transA, *transB, m, n, k, alpha, mhQu, lda, HQ, ldb, beta, innerHess, ldc, len_trans, len_trans);
#endif
#endif
        //diagonalize inner Hessian
        int numeFound;
        int info = diagHessian(innerEigVec, innerEigVal, innerHess, idM, numeFound);
        if( info == 0 ){
            //find max eigenvalue
            //index for sort by absolute and find min/max
            double trce=0.0;
            lambdaMax2 = -2000.0; lambdaMin = 2000.0; 
            for(int i=0;i<idM;i++){
                trce += innerEigVal[i];
                innerEigVal[i] -= spdOffs;	//remove spd correction
                eigIndx[i] = i;
                //report <<debug(1)<<"Eig "<<i<<" is "<<innerEigVal[i]<<endr;
                if(fabs(innerEigVal[i]) > lambdaMax2) lambdaMax2 = fabs(innerEigVal[i]);
                //if(innerEigVal[i] > lambdaMax) lambdaMax = innerEigVal[i];
                if(innerEigVal[i] < lambdaMin) lambdaMin = innerEigVal[i];
            }
            trce -= idM * spdOffs;	//get trace without offsets
            absSort(innerEigVec, innerEigVal, eigIndx, idM);
            lambdaMax = fabs(innerEigVal[_rfM-1]);
            //for(int m=0;m<idM;m++) trce += innerEigVal[m];
            report.precision(5);
            report <<debug(1)<<"[NormalModeDiagonalize::diagonalize] Iteration = " <<iter<< " Eigenvalue upper bound = "<<
                                lambdaMax2<<" Lower bound = "<<lambdaMin<<" Trace = "<<trce<<" _rfM value = "<<
                                    lambdaMax<<" Number eigs found = "<<numeFound<<endr;
        }else{
            report << error << "[NormalModeDiagonalize::diagonalize] Inner diagonalization failed."<<endr;
        }
        //put \bar{Q}_i into mhQu
        getNewEigs(mhQu, innerEigVec, idM);
        //Normalize
        for(int j=0;j<idM;j++){	//NORM
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
            int incxy = 1;
#endif
            n = _3N;
#if defined(HAVE_LAPACK) 
            norm = dnrm2_(&n, &mhQu[j*_3N], &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
            norm = dnrm2_(n, &mhQu[j*_3N], incxy);
#endif
#endif
            if (norm == 0.0) report << error << "[NormalModeDiagonalize] Vector with zero norm."<<endr;
            w[j] = 1.0 / norm;	//re-use w to fix HQ without N^2*m re-calc
            for(int i=0;i<_3N;i++) mhQu[i + j*_3N] /= norm;	//TRY FIXING HQ FOR THE DIAG MMATRIX OF NORMS OR JUST DO \bar{Q}H\bar{Q} (BUT MORE EXPENSIVE)
        }
        //######## Test end and go ##########################################
        if(lambdaMax < traceThreshold) break;
        //######## or find Q_{i+1} ##########################################
        else{
            //removeBondEig();
            traceThreshold = rediagHyst;	//set tighter threshold (hysteresis)
                //##steepest descent method##
                //removeBondEig();	works here!!!!!!!!!!!
                //Update HQ to \bar{Q}H 
                for(int i=0;i<idM;i++)
                    for(int j=0;j<idM;j++) innerEigVec[j + i*idM] *= w[i];
                getNewEigs(HQ, innerEigVec, idM);	
                //Find T1=(I-\bar{Q}_i\bar{Q}_i^T)H\bar{Q}_i, all ops N*m^2
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
#if defined(HAVE_LAPACK)   
                dgemm_ (transA, transB, &m, &n, &k, &alpha, mhQu, &lda, HQ, &ldb, &beta, tempMxM, &ldc);
                transA = "N"; transB = "N"; m = _3N; n = idM; k = idM; lda = _3N; ldb = idM; ldc = _3N; alpha = -1.0; beta = 1.0;
                for(int i=0;i<_3N*idM;i++) T1[i] = HQ[i];	//copy HQ to T1
                dgemm_ (transA, transB, &m, &n, &k, &alpha, mhQu, &lda, tempMxM, &ldb, &beta, T1, &ldc);
                //Find C_i = \bar{Q}H T1
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (transA, transB, &m, &n, &k, &alpha, HQ, &lda, T1, &ldb, &beta, tempMxM, &ldc);
                //copy -C_ii to w_i
                for(int i=0;i<idM;i++) w[i] = -tempMxM[i*idM + i];
                //Find E_i=HQ^T (I-\bar{Q}_i\bar{Q}_i^T)H T1, N^2*m for first op, N*m^2 after
                transA = "N"; transB = "N"; m = _3N; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = _3N; alpha = 1.0; beta = 0.0; //*
                dgemm_ (transA, transB, &m, &n, &k, &alpha, hsn.hessM, &lda, T1, &ldb, &beta, temp3NxM, &ldc);	//(H T1)
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (transA, transB, &m, &n, &k, &alpha, mhQu, &lda, temp3NxM, &ldb, &beta, tempMxM, &ldc);	//\bar{Q}^T (H T1)
                transA = "N"; transB = "N"; m = _3N; n = idM; k = idM; lda = _3N; ldb = idM; ldc = _3N; alpha = -1.0; beta = 1.0;
                //No need to copy, (H T1) already in temp3NxM
                dgemm_ (transA, transB, &m, &n, &k, &alpha, mhQu, &lda, tempMxM, &ldb, &beta, temp3NxM, &ldc);	//(H T1) - \bar{Q}(\bar{Q}^T (H T1))
                //\bar{Q}H ((H T1) - \bar{Q}(\bar{Q}^T (H T1)))
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (transA, transB, &m, &n, &k, &alpha, HQ, &lda, temp3NxM, &ldb, &beta, tempMxM, &ldc);
#else
#if defined(HAVE_SIMTK_LAPACK)
                int len_trans = 1;
                dgemm_ (*transA, *transB, m, n, k, alpha, mhQu, lda, HQ, ldb, beta, tempMxM, ldc, len_trans, len_trans);
                transA = "N"; transB = "N"; m = _3N; n = idM; k = idM; lda = _3N; ldb = idM; ldc = _3N; alpha = -1.0; beta = 1.0;
                for(int i=0;i<_3N*idM;i++) T1[i] = HQ[i];	//copy HQ to T1
                dgemm_ (*transA, *transB, m, n, k, alpha, mhQu, lda, tempMxM, ldb, beta, T1, ldc, len_trans, len_trans);
                //Find C_i = \bar{Q}H T1
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (*transA, *transB, m, n, k, alpha, HQ, lda, T1, ldb, beta, tempMxM, ldc, len_trans, len_trans);
                //copy -C_ii to w_i
                for(int i=0;i<idM;i++) w[i] = -tempMxM[i*idM + i];
                //Find E_i=HQ^T (I-\bar{Q}_i\bar{Q}_i^T)H T1, N^2*m for first op, N*m^2 after
                transA = "N"; transB = "N"; m = _3N; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = _3N; alpha = 1.0; beta = 0.0; //*
                dgemm_ (*transA, *transB, m, n, k, alpha, hsn.hessM, lda, T1, ldb, beta, temp3NxM, ldc, len_trans, len_trans);	//(H T1)
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (*transA, *transB, m, n, k, alpha, mhQu, lda, temp3NxM, ldb, beta, tempMxM, ldc, len_trans, len_trans);	//\bar{Q}^T (H T1)
                transA = "N"; transB = "N"; m = _3N; n = idM; k = idM; lda = _3N; ldb = idM; ldc = _3N; alpha = -1.0; beta = 1.0;
                //No need to copy, (H T1) already in temp3NxM
                dgemm_ (*transA, *transB, m, n, k, alpha, mhQu, lda, tempMxM, ldb, beta, temp3NxM, ldc, len_trans, len_trans);	//(H T1) - \bar{Q}(\bar{Q}^T (H T1))
                //\bar{Q}H ((H T1) - \bar{Q}(\bar{Q}^T (H T1)))
                transA = "T"; transB = "N"; m = idM; n = idM; k = _3N; lda = _3N; ldb = _3N; ldc = idM; alpha = 1.0; beta = 0.0;
                dgemm_ (*transA, *transB, m, n, k, alpha, HQ, lda, temp3NxM, ldb, beta, tempMxM, ldc, len_trans, len_trans);
#endif
#endif
                //make w_i = -C_ii / E_ii
                for(int i=0;i<idM;i++){
                    if(tempMxM[i*idM + i] > 1e-6) w[i] /= tempMxM[i*idM + i]; //MACHINE PRECISION CHECK HERE? ::was e-6
                    else w[i] = 0.0; //vector already an eigenvector if here 
                }
                //Make Q_{i+1} = (\bar{Q}_i + T1 W), normalize? N*m ops
                for(int j=0;j<idM;j++){
                    for(int i=0;i<_3N;i++) mhQu[i + j*_3N] += T1[i + j*_3N] * w[j];
                    //DONT BOTHER WITH NORM HERE AS NEED TO NORM AFTER REDIAG ANYWAY
                }
                //Remove bond eigs (only from \Delta )
                //removeBondEig();
        }
    }
    if(iter >= maxIter){
        report << hint << "[NormalModeDiagonalize] Re-diagonalize exceeded maximum iterations."<<endr;
        iter = 0;
    }
    //
    return iter;
  }

  void NormalModeDiagonalize::removeBondEig(){
    int a1, a2, offs1, offs2;
    Vector3D aDiff;
    Real dotProd;

    for (unsigned int i=0; i<app->topology->bonds.size(); i++){
        //find position vectors and find the normed difference
        a1 = app->topology->bonds[i].atom1; a2 = app->topology->bonds[i].atom2;
        if((app->topology->atoms[a1].name.c_str())[0] == 'H' || (app->topology->atoms[a2].name.c_str())[0] == 'H'){
            //report << hint << "[NormalModeDiagonalize] i= "<<i<<" a1= "<<a1<<" a2 = "<<a2<<" name1= "<<app->topology->atoms[a1].name<<" name2= "<<app->topology->atoms[a2].name<<endr;
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
                //report << hint << "[NormalModeDiagonalize] dp1 = "<<dotProd<<" dp2 = "<<dotProd2<<endr;
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
        if(norm == 0.0) report << error << "[NormalModeDiagonalize] Vector with zero norm."<<endr;
        for(int i=0;i<_3N;i++) mhQu[i + j*_3N] /= norm; 
    }
  }

  //product of eigenvectors and diagonalized 'inner' hessian
  void NormalModeDiagonalize::getNewEigs(double *eigVec, double *innerEigVec, int rfM){
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
    char *transA = "N"; char *transB = "N";
    int m = _3N; int n = rfM; int k = rfM;
    int lda = _3N; int ldb = rfM; int ldc = _3N;
    double alpha = 1.0;	double beta = 0.0;
#endif
    double *tempMat = new double[rfM*_3N];

    //
#if defined(HAVE_LAPACK)
    dgemm_ (transA, transB, &m, &n, &k, &alpha, eigVec, &lda, innerEigVec, &ldb, &beta, tempMat, &ldc);
    //dgemm_ (transA, transB, &m, &n, &k, &alpha, origEigVec, &lda, innerEigVec, &ldb, &beta, eigVec, &ldc);
#else
#if defined(HAVE_SIMTK_LAPACK)
    int len_transa = 1;	int len_transb = 1;
    dgemm_ (*transA, *transB, m, n, k, alpha, eigVec, lda,
                        innerEigVec, ldb, beta, tempMat, ldc, len_transa, len_transb);
#endif
#endif
    //copy back to eigVec
    for(int i=0;i<rfM*_3N;i++) eigVec[i] = tempMat[i];
    delete [] tempMat;
  }

//********************************************************************************************************************************************

  //NVE setps for Hessian average
  void NormalModeDiagonalize::VerletAverage(){
    Real h = avStep * Constant::INV_TIMEFACTOR;

    myAveragePos = *&app->positions;
    numAveragePos = 1;
    calculateForces();
    //Verlet prop.
    for(int nos=1;nos<noAvStep;nos++){
        report <<debug(5)<<"[NormalModeDiagonalize::run] averaging step = "<<nos<<" step size = "<<avStep<<endr;
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

  void NormalModeDiagonalize::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("averageSteps", Value(noAvStep,ConstraintValueType::NotNegative()),1,Text("Hessian averaged over number of steps.")));
    parameters.push_back(Parameter("avStepSize",Value(avStep,ConstraintValueType::NotNegative()),1.0,Text("Step size for Hessian averaging.")));
    parameters.push_back(Parameter("reDiagFrequency", Value(rediagCount,ConstraintValueType::NotNegative()),0,Text("Frequency of re-diagonalization (steps).")));
    parameters.push_back(Parameter("fullDiag",Value(fullDiag,ConstraintValueType::NoConstraints()),false,Text("Full diagonalization?")));
    parameters.push_back(Parameter("removeRand",Value(removeRand,ConstraintValueType::NoConstraints()),false,Text("Remove last random perturbation?")));
    parameters.push_back(Parameter("minSteps", Value(minSteps,ConstraintValueType::NotNegative()),0,Text("Max. number of minimizer steps.")));
    parameters.push_back(Parameter("minLim",Value(minLim,ConstraintValueType::NotNegative()),1.0,Text("Minimization limit kcal mol^{-1}.")));
    parameters.push_back(Parameter("rediagThresh",Value(rediagThresh,ConstraintValueType::NotNegative()),0.0,Text("Re-diagonaliztion threshold.")));
    parameters.push_back(Parameter("rediagHyst",Value(rediagHyst,ConstraintValueType::NotNegative()),0.0,Text("Re-diagonaliztion threshold with hysteresis.")));
    parameters.push_back(Parameter("spdOff",Value(spdOff,ConstraintValueType::NotNegative()),200.0,Text("Factor to ensure Hessisn SPD.")));
    parameters.push_back(Parameter("maxIterations", Value(maxIterations,ConstraintValueType::NotNegative()),1000,Text("Maximum re-diagonalization iterations.")));
    parameters.push_back(Parameter("removeBondedEigs",Value(removeBondedEigs,ConstraintValueType::NoConstraints()),false,Text("Remove bonded eigenvectors?")));
  }

  MTSIntegrator* NormalModeDiagonalize::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeDiagonalize(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],
        values[8],values[9],values[10],values[11],values[12],fg,nextIntegrator);
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
