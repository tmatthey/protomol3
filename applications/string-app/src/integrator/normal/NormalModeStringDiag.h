/*  -*- c++ -*-  */
#ifndef NORMALMODESTRINGDIAG_H
#define NORMALMODESTRINGDIAG_H

#include <protomol/integrator/MTSIntegrator.h>
#include <src/integrator/normal/StringNormalModeUtilities.h>
#include <protomol/integrator/hessian/Hessian.h>

#include <protomol/type/Vector3DBlock.h>

#include <protomol/base/TimerStatistic.h>

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //__________________________________________________ NormalModeStringDiag
  class NormalModeStringDiag :
    public MTSIntegrator, public StringNormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeStringDiag();
    NormalModeStringDiag(int cycles, int avs, Real avss, int redi, bool fDiag,
                          bool rRand, int mins, Real minl, Real redn,
                          Real redhy, Real spd, int maxi, bool rBond, int doM,
                          ForceGroup *overloadedForces,
                          StandardIntegrator *nextIntegrator);
    ~NormalModeStringDiag(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeStringDiag
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void drift();
    void VerletAverage();
    void utilityCalculateForces();
    void removeBondEig();
     
  private:
  public:
     int GetDimension() { return _rfM; }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 14;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);
  protected:
    //virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(const std::vector<Value>& values,
                                  ForceGroup *fg,
                                  StandardIntegrator *nextIntegrator) const;
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class StringNormalModeUtilities
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:  
    //void getNewEigs(double *eigVec, double *innerEigVec, int rfM);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Vector3DBlock diagAt, myAveragePos;
    int numAveragePos;
    Hessian hsn;
    double *eigVal, *innerEigVec, *innerEigVal, *innerHess;
    double *T1, *HQ, *tempMxM, *temp3NxM, *w;
    int *eigIndx;
    double *orig_Eig_Vec;
    bool eigAlloc, firstDiag, fullDiag, removeRand;
    int noAvStep;
    Real avStep;
    int rediagCount, nextRediag;
    int minSteps;
    Real minLim;
    bool validMaxEigv;
    StringNormalModeUtilities *myNextNormalMode, *myLastNormalMode;
    int forceCalc;
    Real lastLambda;
    Real rediagThresh, rediagHyst, spdOff;
    int maxIterations;
    int _idM;
    //Diagnostic data
    Timer rediagTime, hessianTime;
    int hessianCounter, rediagCounter, rediagIters, rediagUpdateCounter;
    //Trace enhancements
    bool removeBondedEigs;

    int doMin;

  public:

    int doIntegrate;
    std::string filename;
    Real dih_value;
  };
}

#endif


