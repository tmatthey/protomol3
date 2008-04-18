/*  -*- c++ -*-  */
#ifndef NORMALMODEDIAGONALIZE_H
#define NORMALMODEDIAGONALIZE_H

#include <protomol/integrator/MTSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/integrator/hessian/Hessian.h>

#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //__________________________________________________ NormalModeDiagonalize
  class NormalModeDiagonalize : public MTSIntegrator, public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeDiagonalize();
    NormalModeDiagonalize(int cycles, int avs, Real avss, int redi, bool fDiag, bool rRand,
                                int mins, Real minl, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    ~NormalModeDiagonalize(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeDiagonalize
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void drift();
    void VerletAverage();
    void utilityCalculateForces();

  private:
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 8;}
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
    virtual MTSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Vector3DBlock diagAt;
    Hessian hsn;
    double *eigVal, *innerEigVec, *innerEigVal, *innerHess, *origEigVec;
    int *eigIndx;
    bool eigAlloc, firstDiag, fullDiag, removeRand;
    int noAvStep;
    Real avStep;
    int rediagCount, nextRediag;
    int minSteps;
    Real minLim;
    bool validMaxEigv;
    NormalModeUtilities *myNextNormalMode, *myLastNormalMode;
    int forceCalc;
    Real lastLambda;

  };

}

#endif


