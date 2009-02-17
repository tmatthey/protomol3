/*  -*- c++ -*-  */
#ifndef OPENMMINTEGRATOR_H
#define OPENMMINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

#if defined (HAVE_OPENMM)

#include "OpenMMContext.h"
#include "Vec3.h"
#include "State.h"
#include "System.h"
#include "HarmonicBondForce.h"
#include "HarmonicAngleForce.h"
#include "LangevinIntegrator.h"
#include "OpenMMContext.h"
#include "NonbondedForce.h"

#endif

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class Vector3DBlock;

  //____ OpenMMIntegrator

  class OpenMMIntegrator : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OpenMMIntegrator();
    OpenMMIntegrator(Real timestep, Real LangevinTemperature,
                              Real gamma, int seed, 
                              bool bond, bool angle, bool nonbonded,  //openMM forces
                              ForceGroup *overloadedForces);
    ~OpenMMIntegrator();

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

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Real myLangevinTemperature;
    Real myGamma;
    int mySeed;
    //

#if defined (HAVE_OPENMM)

    OpenMM::System *system;
    OpenMM::HarmonicBondForce* bonds;
    OpenMM::HarmonicAngleForce* angles;
    OpenMM::NonbondedForce * nonbonded;
    OpenMM::LangevinIntegrator *integrator;
    OpenMM::OpenMMContext *context;
    vector<OpenMM::Vec3> openMMpositions, openMMvelocities, openMMforces;

#endif

    //OpenMM force parameters
    bool HarmonicBondForce, HarmonicAngleForce, NonbondedForce;

  };
  //____ INLINES

  // Constants for GROMACS-PROTOMOL conversion
  namespace Constant {

    extern const Real ANGSTROM_NM;
    extern const Real NM_ANGSTROM;
    extern const Real KCAL_KJ;
    extern const Real KJ_KCAL;
    extern const Real FS_PS;

  }

}
#endif
