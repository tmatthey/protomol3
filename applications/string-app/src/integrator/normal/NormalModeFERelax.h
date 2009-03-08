/*  -*- c++ -*-  */
#ifndef NORMALMODEFERELAX_H
#define NORMALMODEFERELAX_H

#include <protomol/integrator/MTSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //__________________________________________________ NormalModeFERelax
  class NormalModeFERelax : public MTSIntegrator, public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeFERelax();
    NormalModeFERelax(int cycles, Real minimlim, bool rediag, bool simplemin, bool dm, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    ~NormalModeFERelax(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeFERelax
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:  
    void utilityCalculateForces();
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 5;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp* app);
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
    int avItrs, itrs, avMinForceCalc, numSteps;

  private:
    int minCount, forceCalc;
    Real minLim;
    NormalModeUtilities *myPreviousNormalMode;
    Real lastLambda;
    bool reDiag, simpleMin;
    Real randStp;
    bool doMin;

    Vector3DBlock save_pos;

  };

}

#endif


