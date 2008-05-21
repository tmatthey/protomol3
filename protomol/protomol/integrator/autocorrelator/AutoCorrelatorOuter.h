/*  -*- c++ -*-  */
#ifndef AUTOCORRELATOROUTER_H
#define AUTOCORRELATOROUTER_H

#include <protomol/integrator/MTSIntegrator.h>

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;
  class DCDTrajectoryWriter;

  class AutoCorrelatorOuter : public MTSIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    AutoCorrelatorOuter();
    AutoCorrelatorOuter(int cycles, string dcdfile, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    ~AutoCorrelatorOuter(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 2;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class MTSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    DCDTrajectoryWriter* myWriter;
    string myDCDFile;
    int numNonWaters;
    Vector3DBlock* nonWater;

  };


}



#endif
