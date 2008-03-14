/*  -*- c++ -*-  */
#ifndef OUTPUTCOLLECTION_H
#define OUTPUTCOLLECTION_H

#include <list>

namespace ProtoMol {
  class Output;
  class OutputFactory;
  class ProtoMolApp;
  //____ OutputCollection

  /**
     Container class for Output objects invoked at application level.
   */
  class OutputCollection  {
    friend class OutputFactory;

    typedef std::list<Output *> Container;
    typedef std::list<Output *>::iterator iterator;

  public:
    typedef std::list<Output *>::const_iterator const_iterator;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputCollection();
    ~OutputCollection();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputCollection
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const ProtoMolApp *app);

    ///< Initialize all Output objects
    void run(int step);

    ///< Invoke all Output objects with run()
    void updateNextStep(int step);
    void finalize(int step);

    ///< Finalize all Outout objects
    int getNext() const;

    ///< Add new Output object to the collection
    void adoptOutput(Output *output);

    /// Iterators, const
    const_iterator begin() const {return myOutputList.begin();}
    const_iterator end()   const {return myOutputList.end();}

  private:
    /// Iterators
    iterator begin()       {return myOutputList.begin();}
    iterator end()         {return myOutputList.end();}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Container myOutputList;

    const ProtoMolApp *app;
  };
}
#endif
