/*  -*- c++ -*-  */
#ifndef OPENMMINTEGRATOR_H
#define OPENMMINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

#if defined (HAVE_OPENMM)

#include "OpenMMContext.h"
#include "Vec3.h"
#include "State.h"
#include "System.h"

//forces
#include "HarmonicBondForce.h"
#include "HarmonicAngleForce.h"
#include "NonbondedForce.h"
#include "RBTorsionForce.h"
#include "PeriodicTorsionForce.h"
#include "GBSAOBCForce.h"
#include "CMMotionRemover.h"

//integrators
#include "Integrator.h"
#include "LangevinIntegrator.h"
//other
#include "OpenMMContext.h"

#endif

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class Vector3DBlock;

  //____ OpenMMIntegrator

  class OpenMMIntegrator : public STSIntegrator {

  private:
    static const Real FudgeQQ;
    static const Real FudgeLJ;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OpenMMIntegrator();
    OpenMMIntegrator(Real timestep, ForceGroup *overloadedForces);
    ~OpenMMIntegrator();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;
    virtual unsigned int getParameterSize() const{return 12;}

  protected:
    virtual void setupValues(std::vector<Value> &params);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OpenMMIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void getObcScaleFactors(vector<Real>& scaleFactors);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  protected:
    Real myLangevinTemperature;
    Real myGamma;
    int mySeed;
    //

#if defined (HAVE_OPENMM)

    OpenMM::System *system;
    OpenMM::HarmonicBondForce* bonds;
    OpenMM::HarmonicAngleForce* angles;
    OpenMM::RBTorsionForce* RBDihedral;
    OpenMM::PeriodicTorsionForce* PTorsion;
    OpenMM::NonbondedForce * nonbonded;
    OpenMM::GBSAOBCForce* gbsa;
    OpenMM::Integrator *integrator;//Langevin
    OpenMM::OpenMMContext *context;
    vector<OpenMM::Vec3> openMMpositions, openMMvelocities, openMMforces;

#endif

    //OpenMM force/integrator parameters
    bool HarmonicBondForce, HarmonicAngleForce, NonbondedForce, RBDihedralForce, PeriodicTorsion;
    int  myIntegratorType;
    Real myGBSAEpsilon, myGBSASolvent; 

  };
  //____ INLINES


}
#endif
