/*  -*- c++ -*-  */
#ifndef OUTPUT_FAH_GUI_H
#define OUTPUT_FAH_GUI_H
#ifdef HAVE_LIBGUI

#include <protomol/output/Output.h>
#include <string>

namespace FAH {
  class GUIServer;
}

namespace ProtoMol {
  class Configuration;

  //____ OutputFAHGUI
  class OutputFAHGUI : public Output {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputFAHGUI();
    OutputFAHGUI(const std::string &name, int freq);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Output *doMake(const std::vector<Value> &values) const;
    virtual void doInitialize();
    virtual void doRun(int step);
    virtual void doFinalize(int);
    virtual bool isIdDefined(const Configuration *config) const;
    virtual bool addDoKeyword() const {return false;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &) const;
    virtual bool adjustWithDefaultParameters(std::vector<Value> &values,
                                             const Configuration *config) const;

  private:
    void setCoords();
    void setBonds();
    void setAtoms();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    std::string name;
    FAH::GUIServer *server;
  };
}

#endif // HAVE_LIBGUI
#endif
