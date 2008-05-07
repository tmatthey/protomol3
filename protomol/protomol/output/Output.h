/*  -*- c++ -*-  */
#ifndef OUTPUT_H
#define OUTPUT_H

#include <protomol/base/Makeable.h>

namespace ProtoMol {
  class ProtoMolApp;

  //____ Output

  /**
     Base class of all Output classes to dump data at a given
     frequency.  The actual output frequency is defined by global
     output frequency times the the frequency of the concrete
     class. The global output frequency is also  used to define the
     number of steps to run the integrator.  It keeps pointers to the
     topology, positions, velocities, and integrator.  Output objects
     are aggregated in OutputCollection, which  invokes them one by
     one at application level.  If you only need to print some output
     to a file you should rather inherit from OutputFile, then inherit
     directly from Output.
   */
  class Output : public Makeable<Output> {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Output();
    Output(int freq);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const ProtoMolApp *app);

    ///< To initialize the object, before the simulation starts.

    void run(int step);

    //< Called at each step (e.g., printing total energy on the screen),   
    //< takes care of the output frequency.

    void finalize(int step);

    //< At the end of the simulation (e.g., writing final positions), and
    //< calls first run() to ensure that run is called for the last
    //< step, if needed.

    Output *make(const std::vector<Value> &values) const;

    ///< Factory method to create a complete output object from its prototype

    virtual bool isIdDefined(const Configuration *config) const;

    ///< Should return true if the concrete object is defined/specified in
    ///< Configuration by the user. Normally if gedId() has a valid value
    ///< in Configuration.

    virtual bool addDoKeyword() const {return true;}

    ///< Defines if the output object supports do<getId()> to enable or disable
    ///< the output.

    int getFirstStep() const {return myFirstStep;}
    int getLastStep() const {return myLastStep;}
    int getOutputFreq() const {return myOutputFreq;}
    int getNext() const {return myNextStep;}
    bool first() const {return myFirst;}
    void updateNextStep(int step);

  private:
    virtual void doInitialize() = 0;

    ///< Hook method of initialize, implemented in the concrete class
    virtual void doRun(int step) = 0;

    ///< Hook method of run, implemented in the concrete class
    virtual void doFinalize(int step) = 0;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getScope() const {return scope;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string scope;
    const ProtoMolApp *app;

  private:
    int myFirstStep;
    int myNextStep;
    int myLastStep;
    bool myFirst;

  protected:
    int myOutputFreq;       ///< Output freqeuncy

  };
}
#endif
