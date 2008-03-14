/*  -*- c++ -*-  */
#ifndef NORMALMODEBROWNIAN_H
#define NORMALMODEBROWNIAN_H

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

//####diagnostics
#include <protomol/io/XYZTrajectoryWriter.h>

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;

  //____ NormalModeBrownian
  class NormalModeBrownian : public STSIntegrator, public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeBrownian();
    NormalModeBrownian(
      Real timestep, int firstmode, int nummode, Real gamma, int seed,
      Real temperature,
      //####added avff, inff for diagnostics
      std::string avff, std::string inff,
      ForceGroup *overloadedForces);
    ~NormalModeBrownian();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeBrownian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
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
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeBrownian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void forceProjection();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    Real randStp;
    int aveForceCount;
    //####diagnostics
    // mean force output:
    std::string avForceFile;
    // instantaneous force output:
    std::string inForceFile;

    XYZTrajectoryWriter *myWriter;
    XYZTrajectoryWriter *myWriter2;
  };
}

#endif

