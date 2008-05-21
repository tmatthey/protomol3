/*  -*- c++ -*-  */
#ifndef AUTOCORRELATORINNER_H
#define AUTOCORRELATORINNER_H


#include <string>
#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>


namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class DCDTrajectoryWriter;

  class AutoCorrelatorInner : public LeapfrogIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    AutoCorrelatorInner();
    AutoCorrelatorInner(Real timestep, string dcdfile, ForceGroup *overloadedForces);
    ~AutoCorrelatorInner();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 2;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);

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
    int getNumNonWaters() {return numNonWaters;}

  private:
    void writeDCD();
    DCDTrajectoryWriter* myWriter;
    string myDCDFile;
    int numNonWaters;
    Vector3DBlock* nonWater;
  };

}


#endif
