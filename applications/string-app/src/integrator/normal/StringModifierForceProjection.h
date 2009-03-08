/*  -*- c++ -*-  */
#ifndef MODIFIERFORCEPROJECTION_H
#define MODIFIERFORCEPROJECTION_H

#include <protomol/modifier/Modifier.h>
#include <src/integrator/normal/StringNormalModeUtilities.h>

namespace ProtoMol {
  //____ StringModifierForceProjection
  class StringModifierForceProjection : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    StringModifierForceProjection(StringNormalModeUtilities *i) :
      Modifier(Constant::MAX_INT - 400), myTheIntegrator(i) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
    virtual std::string getIdNoAlias() const {
      return std::string("ForceProjection");
    };

  private:
    virtual void doExecute(Integrator *i) {
      myTheIntegrator->forceProjection();
    }
    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    StringNormalModeUtilities *myTheIntegrator;
  };
}
#endif /* MODIFIER_H */
