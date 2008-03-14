/*  -*- c++ -*-  */
#ifndef NORMALMODEMINIMIZER_H
#define NORMALMODEMINIMIZER_H

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;

  //____ NormalModeMinimizer
  class NormalModeMinimizer : public STSIntegrator,
    public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeMinimizer();
    NormalModeMinimizer(Real timestep, int firstmode, int nummode, Real gamma,
                        int seed, Real temperature,
                        Real minimlim, bool rforce, bool rediag, bool simplemin,
                        ForceGroup *overloadedForces);
    ~NormalModeMinimizer();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeMinimizer
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void utilityCalculateForces();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;

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
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    int avItrs, itrs, avMinForceCalc, numSteps;

  private:
    int minCount, forceCalc;
    Real minLim;
    bool randforce;
    NormalModeUtilities *myPreviousNormalMode;
    Real lastLambda;
    bool reDiag, simpleMin;
    Real randStp;
  };
}

#endif

