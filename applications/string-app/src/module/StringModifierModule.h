#ifndef STRINGMODIFIER_MODULE_H
#define STRINGMODIFIER_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>

#include <string>

namespace ProtoMol {
  class ProtoMolApp;

  declareInputValue(InputRemoveAngularMomentum, INT, NOCONSTRAINTS);
  declareInputValue(InputRemoveLinearMomentum, INT, NOCONSTRAINTS);
  declareInputValue(InputRattle, BOOL, NOCONSTRAINTS);
  declareInputValue(InputRattleEpsilon, REAL, NOTNEGATIVE);
  declareInputValue(InputRattleMaxIter, INT, NOTNEGATIVE);
  declareInputValue(InputShake, BOOL, NOCONSTRAINTS);
  declareInputValue(InputShakeEpsilon, REAL, NOTNEGATIVE);
  declareInputValue(InputShakeMaxIter, INT, NOTNEGATIVE);
  declareInputValue(InputShadow, BOOL, NOCONSTRAINTS);
  declareInputValue(InputShadowOrder, INT, NOTNEGATIVE);
  declareInputValue(InputShadowFreq, INT, NOTNEGATIVE);
  declareInputValue(InputHarmNM, INT, NOTNEGATIVE);
  declareInputValue(InputOutputIntermediate, INT, NOTNEGATIVE);
  declareInputValue(InputOutputIntermediateStart, INT, NOTNEGATIVE);

  class StringModifierModule : public Module {
  public:
    const std::string getName() const {return "Modifier";}
    int getPriority() const {return 0;}
    const std::string getHelp() const {return "";}

    void init(ProtoMolApp *app);
    void postBuild(ProtoMolApp *app);
    void addModifiers(ProtoMolApp *app);
  };
};

#endif // STRINGMODIFIER_MODULE_H
