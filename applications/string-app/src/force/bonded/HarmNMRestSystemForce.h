/* -*- c++ -*- */
#ifndef HARMNMRESTSYSTEMFORCE_H
#define HARMNMRESTSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/io/EigenvectorReader.h>
#include <protomol/io/PDBReader.h>
#include <src/base/CommonSpace.h>

#include <string>

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#else
#if defined (HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

#define NUMPOINTS 2
#define THRES_HOLD 0.95
#define INC_KAPPA 0.5

using namespace ProtoMol::Report;
using std::cout;

namespace ProtoMol {

  //____ HarmNMRestSystemForce

  template<class TBoundaryConditions>
  class HarmNMRestSystemForce :
       public SystemForce, public NormalModeUtilities {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

      HarmNMRestSystemForce() : SystemForce(), NormalModeUtilities() {

         kappaScale=1.0;
         doReadParams = 0;
         firstMode = 1;
         numMode = 1;
         eigvecfile = "";
         eigvecfileB = "";
         pdbFileA = "";
         pdbFileB = "";
         target_point = 0;
         target_steps = 0;
         doCommonSpace = 0;
         target_vector = 0;
         massMatrix = 0;
         eigVals = NULL;
         tC = NULL;
         mdiff = NULL; /* array for holding the diff between modes */
         refPos = new Vector3DBlock;
         oldRefPos = new Vector3DBlock;
         refPosB = new Vector3DBlock;

         tmpVec = new Vector3DBlock;

         cHat = NULL;
         cP = NULL;

         step = 0;

         mhQu = NULL;

         reparam_target = 0;

      }

      HarmNMRestSystemForce(Real ks, int fm, int nm, int doR, std::string evfile, std::string evfile_B, std::string pdbA, std::string pdbB, int t_p, int t_s, int do_c, Real n_s, int t_v, Real max_k, int start_common, Real th_hold, int rp_target) : 
               SystemForce(), NormalModeUtilities(fm, nm,91.0,1234,300),kappaScale(ks), doReadParams(doR),
               eigvecfile(evfile), eigvecfileB(evfile_B), pdbFileA(pdbA),pdbFileB(pdbB), target_point(t_p), target_steps(t_s),
               doCommonSpace(do_c), norm_epsilon(n_s), target_vector(t_v), max_kappa(max_k), start_common_mode(start_common),
               thres_hold(th_hold), reparam_target(rp_target) {

         massMatrix = 0;
         eigVals = NULL;
         tC = NULL;
         refPos = new Vector3DBlock;
         oldRefPos = new Vector3DBlock;
         refPosB = new Vector3DBlock;
         tmpVec = new Vector3DBlock;

         cHat = NULL;
         cP = NULL;
         mdiff = NULL; /* array for holding the diff between modes */

         step = 0;
         mhQu = NULL;

      }
      
     virtual ~HarmNMRestSystemForce() {
         if (refPos != NULL) delete refPos;
         if (oldRefPos != NULL) delete oldRefPos;
         if (refPosB != NULL) delete refPosB;
         if (tmpVec != NULL) delete tmpVec;
         if (tC != NULL) delete [] tC;

         if (cHat != NULL) delete [] cHat;
         if (cP != NULL) delete [] cP;
         if (mdiff != NULL) delete [] mdiff;

     }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New method of HarmNMRestSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    void initializeDataStructures(const Vector3DBlock *positions,
          const GenericTopology *topo, Vector3DBlock *forces, int nm_flags) ;
    Vector3DBlock *subspaceForceCartProj(double *tempC, Vector3DBlock *iPforce);
    Vector3DBlock *subspacePositionHalfProj(double *tempC, Vector3DBlock *iPforce);
    Vector3DBlock *subspaceKappaForce(Vector3DBlock *force, Vector3DBlock *iPforce);

    void SetEigenvalues(); //function for copying the eigenvectors from integrator

    void SetCHat(double *ch, int doZero);

    void SetStep(int s_t) { step = s_t; }

    void Check_Kappa_And_Upgrade();

    void Set_Next_Target() {
           cout<<" HarmNMForce : updating cHat..._rfM "<<_rfM<<endl;
           for(int i=0;i<(_rfM-(firstMode-1));i++) cHat[i] += cP[i];
           for(int i=0;i<(_rfM-(firstMode-1));i++) cout<<"Mode "<<(i+firstMode)<<" cHat target "<<cHat[i]<<endl;
           cout<<endl;
    }
       

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions, Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "HarmNMRestForce";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos){return 0;}

   private:
    virtual Force *doMake(const std::vector<Value> &values) const {
       return new HarmNMRestSystemForce(values[0],values[1],values[2],values[3],values[4],values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16]);
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &parameters) const {
       parameters.push_back
          (Parameter("-scale", Value(kappaScale, ConstraintValueType::NotNegative()), 
                   Text("potential bias constant")));
       parameters.push_back(Parameter("-firstMode", Value(firstMode, ConstraintValueType::NotNegative()),
                   Text("First Mode")));
       parameters.push_back(Parameter("-numberModes", Value(numMode, ConstraintValueType::NotNegative()),
                   Text("No of Modes")));
       parameters.push_back(Parameter("-doReadParams", Value(doReadParams, ConstraintValueType::NotNegative()),0,
                   Text("do read other parameters?"))); 
       parameters.push_back(Parameter("-eigvecfile1", Value(eigvecfile, ConstraintValueType::NoConstraints()),"",
                   Text("eigenvector file")));
       parameters.push_back(Parameter("-eigvecfile2", Value(eigvecfileB, ConstraintValueType::NoConstraints()),"",
                   Text("eigenvector file in state B")));
       parameters.push_back(Parameter("-pdbfileA", Value(pdbFileA, ConstraintValueType::NoConstraints()), "", Text("pdb file A")));
       parameters.push_back(Parameter("-pdbfileB", Value(pdbFileB, ConstraintValueType::NoConstraints()), "", Text("pdb file B")));
       parameters.push_back(Parameter("-point", Value(target_point, ConstraintValueType::NotNegative()), 1, Text("Target point")));
       parameters.push_back(Parameter("-steps", Value(target_steps, ConstraintValueType::NotNegative()), 0, Text("Target point")));
       parameters.push_back(Parameter("-doCommonSpace", Value(doCommonSpace, ConstraintValueType::NoConstraints()), 0, Text("do common space?")));
       parameters.push_back(Parameter("-norm_epsilon", Value(norm_epsilon, ConstraintValueType::NoConstraints()), 0.1, Text("Epsilon for finding closest point")));
       parameters.push_back(Parameter("-target_vector", Value(target_vector, ConstraintValueType::NoConstraints()), 0, Text("target_vector for studying restraint")));
       parameters.push_back(Parameter("-max_kappa", Value(max_kappa, ConstraintValueType::NoConstraints()), 25.0, Text("Maximum kappa for any subspace coordinate")));
       parameters.push_back(Parameter("-start_common", Value(start_common_mode, ConstraintValueType::NoConstraints()), 0.0, Text("If common space, start restraint from this vector")));
       parameters.push_back(Parameter("-threshold_commonspace", Value(thres_hold, ConstraintValueType::NoConstraints()),0.9, Text("Threshold for finding common space")));
       parameters.push_back(Parameter("-reparam_target", Value(reparam_target, ConstraintValueType::NoConstraints()),0, Text("Reparameterization target")));

    }

   private:
     virtual void doSetParameters(std::string &, std::vector<Value> values) {
      kappaScale = values[0];
      firstMode = values[1];
      numMode = values[2];
      doReadParams = values[3];
      eigvecfile = std::string(values[4]);
      eigvecfileB = std::string(values[5]);
      pdbFileA = std::string(values[6]);
      pdbFileB = std::string(values[7]);
      target_point = values[8];
      target_steps = values[9];
      doCommonSpace = values[10];
      norm_epsilon = values[11];
      target_vector = values[12];
      max_kappa = values[13];
      start_common_mode = values[14];
      thres_hold = values[15];
     }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    public:

       Real kappaScale;

       int doReadParams;
 
       std::string eigvecfile, eigvecfileB;

       std::string pdbFileA, pdbFileB;

       int target_point, target_steps;

       int doCommonSpace;

       Real norm_epsilon;

       int target_vector;

       Real max_kappa;

       int start_common_mode;

       Real thres_hold;

       double *eigVals;

       Vector3DBlock *refPos, *oldRefPos, *refPosB;

       Vector3DBlock *tmpVec;


    public:
       // flag for running initialDataStructures first time evaluate is called
       int massMatrix;


       double *tC;

       double *cHat; //holds the target C-values

       double *cP;

       double *mdiff;

       PDBReader pdbReader;

       PDB pdbA;

       PDB pdbB;

       EigenvectorReader eigvecReader;

       EigenvectorInfo eigenvecInfo, eigenvecInfo_B;

       Vector3DBlock zeroV;

       int step;

       int satisfy_maxnorm;

       int next_target_steps;

       int reparam_target;

  };

    template<class TBoundaryConditions>
  inline void HarmNMRestSystemForce<TBoundaryConditions>::initializeDataStructures(const Vector3DBlock *positions,
          const GenericTopology *topo, Vector3DBlock *forces, int nm_flags) {

          //NormalModeUtilities::initialize((int)positions->size(), topo, forces, false);
#if defined(HAVE_LAPACK)
#else
#if defined(HAVE_SIMTK_LAPACK)
#else
    report << error << "Normal Mode integrators require Lapack libraries."<<endr;
#endif
#endif
    //local force pointer
    myForcesP = forces;
    //type
    if(nm_flags & COMPLIMENT_FORCES) complimentForces = true;
    else  complimentForces = false;
    if(nm_flags & GEN_COMP_NOISE) genCompNoise = true;
    else  genCompNoise = false;
    //
    _N = (int)positions->size();
    _3N = 3*_N; //used everywhere, must be initialized

    //setup for auto re-diag by allowing full set of eigenvectors
    if(!numEigvectsu) numEigvectsu = _3N;
    std::cout<<"_N "<<_N<<"_3N "<<_3N<<"numeigvectsu "<<numEigvectsu<<std::endl;
    //first mode?
    if(firstMode < 1 || firstMode > _3N) report << error << "firstmode = "<<firstMode<<", must be > 0 and less than 3N = "<<_3N<<"."<<endr;
    //check numMode
    if(numMode < 1 || numMode > _3N - (firstMode-1)) report << error << "numbermodes = "<<numMode<<", must be > 0 and less than or equal to 3N-(firstmode-1) = "<<_3N-(firstMode-1)<<"."<<endr;
    //
    if(numMode > (int)numEigvectsu - (firstMode-1))  report << error << "Insufficient eigenvectors: numbermodes = "<<numMode<<", must be less than than or equal to m-(firstmode-1) = "<<numEigvectsu-(firstMode-1)<<"."<<endr;
    //number of low frequency modes
    _rfM = numMode + (firstMode-1);
   //temporary mode space variable for intermediate calculations
    tmpC = new double[_3N];
    //define temporary position/force vector
    //tmpFX = new double[_3N];
    //setup sqrt masses and their inverse
    invSqrtMass = new double[_N];
    sqrtMass = new double[_N];
    for(int i=0;i<_N;i++){
        //cout<<"Mass of atom "<<i<<" "<<topo->atoms[i].scaledMass<<endl;
        sqrtMass[i] = sqrt( topo->atoms[i].scaledMass );
        if(sqrtMass[i]) invSqrtMass[i] = 1.0 / sqrtMass[i];
    }


     tC = new double[_3N];


          

    }
          
  template<class TBoundaryConditions>
  inline void HarmNMRestSystemForce<TBoundaryConditions>::evaluate(
      const GenericTopology *topo, const Vector3DBlock *positions,
      Vector3DBlock *forces, ScalarStructure *energies) {


        if (!massMatrix) {

           if (doReadParams) {
                //read input files
                
                //read the pdb files
                if (!pdbReader.open(pdbFileA)) report << error <<"HarmNMRestSystemForce : Cant open PDB file A "<<endr;
                if (!(pdbReader >> pdbA)) report << error <<"HarmNMRestSystemForce : Cant read PDB file A "<<endr;
                pdbReader.close(); 

                if (!pdbReader.open(pdbFileB)) report << error <<"HarmNMRestSystemForce : Cant open PDB file A "<<endr;
                if (!(pdbReader >> pdbB)) report << error <<"HarmNMRestSystemForce : Cant read PDB file A "<<endr;
                pdbReader.close(); 

                refPos->intoAssign(pdbA.coords);
                refPosB->intoAssign(pdbB.coords);

                //read eigenvector files
                if (!eigvecReader.open(eigvecfile)) report << error <<"HarmNMRestSystemForce : Cant open eigenvector file"<<endr;
                if (!(eigvecReader >> eigenvecInfo)) report << error <<"HarmNMRestSystemForce : Cant read eigenvector file"<<endr;
                eigvecReader.close();

                //check if eigvecfileB is available
                if (!eigvecfileB.empty()) {
                   cout<<"Second eigenvector file name "<<eigvecfileB<<endl;
                   if (!eigvecReader.open(eigvecfileB)) report << error <<"HarmNMRestSystemForce : Cant open eigenvector file B"<<endr;
                   if (!(eigvecReader >> eigenvecInfo_B)) report << error <<"HarmNMRestSystemForce : Cant read eigenvector file B"<<endr;
                   eigvecReader.close();
                }

                mhQu = eigenvecInfo.myEigenvectors;
                Q = &mhQu;
                
           }

           if (refPos->empty()) report << error <<"x_0 is not set in  HarmNMRestSystemForce"<<endr;

           initializeDataStructures(positions, topo, forces, false);


           //find C positions
           if (doReadParams) {
               if (doCommonSpace) {
                   //find the common space and assign to Q
                   CommonSpace *myCommonSpace = new CommonSpace();
                   //myCommonSpace->initialize(NUMPOINTS, _3N, _rfM, THRES_HOLD, eigenvecInfo.myEigenvectors);
                   myCommonSpace->initialize(NUMPOINTS, _3N, _rfM, thres_hold, eigenvecInfo.myEigenvectors);
                   vector<double *> eigvectors;
                   double *ev_A = eigenvecInfo.myEigenvectors;
                   double *ev_B = eigenvecInfo_B.myEigenvectors;
                   eigvectors.push_back(ev_A);
                   eigvectors.push_back(ev_B);
                   myCommonSpace->run(eigvectors);
                   if (mhQu == NULL) mhQu = new double[_3N*_3N];
                   double *ev = myCommonSpace->commonBasis;
                   _rfM = myCommonSpace->basis_size;
                   for(int k=0;k<_rfM*_3N;k++) mhQu[k] = ev[k];
                   //if firstmode not set to 1, we have a problem

                   delete(myCommonSpace);

                   if (target_vector > 0) {
                      firstMode = target_vector;
                      _rfM = firstMode;
                   }
                    
                   if (start_common_mode > 0) firstMode = start_common_mode;
               }
                   

              if (cHat == NULL) cHat = new double[_rfM];
              if (cP == NULL) cP = new double[_rfM];

              modeSpaceProj(cP,refPosB,refPos);
              for(int i=0;i<(_rfM-(firstMode-1));i++) {
                   cP[i] /= target_point;
                   //cout<<"target mode "<<cP[i]<<endl;
              }
              for(int i=0;i<(_rfM-(firstMode-1));i++) cHat[i] = 0;
               
              
              if (eigVals == NULL) eigVals = new double[_rfM];
              for(int j=0;j<(_rfM-(firstMode-1));j++) eigVals[j] = kappaScale;

           }
           
           if (refPos->empty()) report << error <<"x_0 is not set in  HarmNMRestSystemForce"<<endr;
           if (eigVals == NULL) report << error <<"Spring constants not set"<<endr;

           if (mdiff == NULL) mdiff = new Real[_rfM];
           massMatrix = 1;
           
           //next_target_steps : the no of steps on which next update in C value will take place.
           next_target_steps = target_steps;


        }else {

        if (!reparam_target) {
        if ((step%next_target_steps) == 0) {
           //update cHat
           Set_Next_Target();
        }

        }

        }

        tmpVec->intoAssign(*positions);

        //mode space projection. In this step, I am doing Q^{T}M^{0.5}(x_c - x_0).
        modeSpaceProj(tC, tmpVec, refPos);

        Check_Kappa_And_Upgrade();

        /*for(int i=0;i<(_rfM-(firstMode-1));i++) {
           cout<<"Mode "<<(i+(firstMode))<<" C-val "<<tC[i]<<endl;
        }*/

       for(int i=0;i<(_rfM-(firstMode-1));i++) mdiff[i] = tC[i] - cHat[i];

       satisfy_maxnorm = 1;
       for(int j=0;j<(_rfM-(firstMode-1));j++) {
           //cout<<"j "<<j<<" "<<mdiff[j]<<endl;
           if (fabs(mdiff[j]) > norm_epsilon) {
             satisfy_maxnorm = 0;
             break; 
           }/* else continue to see if all diffs are below the norm_epsilon */
       }
       if (satisfy_maxnorm) report << plain <<"Max norm has been reached at step "<<step<<endl;

        for(int i=0;i<(_rfM-(firstMode-1));i++) {
           tC[i] = eigVals[i]*(tC[i] - cHat[i]);

        }

        //subspaceForceCartProj(tC,tmpVec); //projection of mean force in x space.
        subspacePositionHalfProj(tC,tmpVec);


        Real harm_pot = 0;

        for(int i=0;i<_3N;i++) harm_pot += tmpVec->c[i]*tmpVec->c[i];

        harm_pot *= 0.5;

	cout<<"HarmNMRestForce : Potential due to harmonic restraint "<<harm_pot<<endl;

        (*energies)[ScalarStructure::OTHER] = harm_pot;

        subspaceKappaForce(tmpVec, tmpVec);
        
        //subspaceForce(tmpVec,tmpVec);
        for(int i=0;i<_3N;i++) forces->c[i] -= tmpVec->c[i];

   }

       template<class TBoundaryConditions>
  inline Vector3DBlock *HarmNMRestSystemForce<TBoundaryConditions>::subspacePositionHalfProj(
           double *tempC, Vector3DBlock *iPos) {

    //The code here is a variant of cartSpaceProj in NormalModeUtilities
        char *transA = "N";
        int m = _3N; int n = _rfM-(firstMode-1); int incxy = 1;
        double alpha = 1.0; double beta = 0.0;

#if defined(HAVE_LAPACK)
    dgemv_ (transA, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, tempC, &incxy, &beta, iPos->c, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
    int len_transa = 1;                         //length of transA
    dgemv_ (*transA, m, n, alpha, &((*Q)[_3N*(firstMode-1)]), m, tempC, incxy, beta, iPos->c, incxy, len_transa);
#endif
#endif
    //
        for( int i=0; i < _3N; i++)
            iPos->c[i] /= sqrtMass[i/3];

    return iPos;

   }

  //Find forces acting inside subspace
  template<class TBoundaryConditions>
  inline Vector3DBlock *HarmNMRestSystemForce<TBoundaryConditions>:: subspaceKappaForce(Vector3DBlock * force, Vector3DBlock * iPforce){
    //
    if (force != iPforce)
      iPforce->intoAssign(*force);
    //calculate M^{1/2}QQ^TM^{-1/2}f using BLAS
    //f'=M^{-1/2}*f
    for( int i=0; i < _3N; i++) {
            iPforce->c[i] *= invSqrtMass[i/3];
    }
    //c=Q^T*M^{-1/2}*f
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
    char *transA = "T";                                                 // Transpose, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM-(firstMode-1); int incxy = 1;     //sizes
#endif
    double alpha = 1.0; double beta = 0.0;
    //
#if defined(HAVE_LAPACK)
    dgemv_ (transA, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, iPforce->c, &incxy, &beta, tmpC, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
    int len_transa = 1;
    dgemv_ (*transA, m, n, alpha, &((*Q)[_3N*(firstMode-1)]), m, iPforce->c, incxy, beta, tmpC, incxy, len_transa);
#endif
#endif

    //Add \kappa term
    for(int j=0;j<_rfM-(firstMode-1);j++) tmpC[j] *= eigVals[j];

    //
    //f''=Qc
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
    char *transB = "N"; /* LAPACK checks only first character N/V */
#endif
    alpha = 1.0;        beta = 0.0;
    //
#if defined(HAVE_LAPACK)
    dgemv_ (transB, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, tmpC, &incxy, &beta, iPforce->c, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
    dgemv_ (*transB, m, n, alpha, &((*Q)[_3N*(firstMode-1)]), m, tmpC, incxy, beta, iPforce->c, incxy, len_transa);
#endif
#endif
    //f'''=M^{1/2}*f''
    for( int i=0; i < _3N; i++) {
            iPforce->c[i] *= sqrtMass[i/3];
    }
    //put back into vector3DBlocks
    //vectTOvector3DBlock(tmpFX, iPforce);
    //delete temporary array
    return iPforce;
  }


      template<class TBoundaryConditions>
   inline  Vector3DBlock *HarmNMRestSystemForce<TBoundaryConditions>::subspaceForceCartProj(
       double *tempC, Vector3DBlock *iPforce) {

       char *transB = "N";
       int m = _3N; int n = _rfM; int incxy = 1;
       double alpha = 1.0; double beta = 0.0;

#if defined(HAVE_LAPACK)
    dgemv_ (transB, &m, &n, &alpha, &((*Q)[_3N*(firstMode-1)]), &m, tempC, &incxy, &beta, iPforce->c, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
    dgemv_ (*transB, m, n, alpha, &((*Q)[_3N*(firstMode-1)]), m, tempC, incxy, beta, iPforce->c, incxy, len_transa);
#endif
#endif
    for( int i=0; i < _3N; i++) {
            iPforce->c[i] *= sqrtMass[i/3];
    }

   return iPforce;


   }


 template<class TBoundaryConditions>
  inline void HarmNMRestSystemForce<TBoundaryConditions>::SetEigenvalues() {


      cout<<"HarmNMForce , CopyEigenvectors _rfM "<<_rfM<<endl;

      if (eigVals == NULL) eigVals = new double[_rfM];

      //for(int j=0;j<(_rfM-(firstMode-1));j++) eigVals[j] = 1.0;
      for(int j=0;j<(_rfM-(firstMode-1));j++) eigVals[j] = kappaScale;

   } 

  template<class TBoundaryConditions>
    inline void HarmNMRestSystemForce<TBoundaryConditions>::SetCHat(double *ch,int doZero) {

       if (cHat == NULL) {
          cHat = new double[_rfM];
       }
       if (!doZero) {
          for (int j=0;j<(_rfM-(firstMode-1));j++) cHat[j] = 0;
       }else {
          for(int j=0;j<(_rfM-(firstMode-1));j++) cHat[j] += ch[j];
       }
   }


   template<class TBoundaryConditions>
   inline void HarmNMRestSystemForce<TBoundaryConditions>::Check_Kappa_And_Upgrade() {
        if (step % 1000 != 0)  return;

        Real mdiff;
        
        for (int j=0 ; j<(_rfM - (firstMode-1));j++) {
           mdiff = fabs(tC[j]-cHat[j]);
           if (fabs(mdiff) > norm_epsilon) {
              if (eigVals[j] < max_kappa) {
                 cout <<"At step "<<step<<" ,incremeting kappa for mode "<<j+firstMode<<" from "<<eigVals[j]<<" to ";
                 eigVals[j] += INC_KAPPA;
                 cout<<eigVals[j]<<endl;
              }
           }
        }
   }
           

}        

#endif
