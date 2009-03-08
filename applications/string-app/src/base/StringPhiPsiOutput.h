/* -*- c++ -*- */
#ifndef STRINGPHIPSIOUTPUT_H
#define STRINGPHIPSIOUTPUT_H

#include <src/base/MyOutput.h>

class StringPhiPsiOutput : public MyOutput {

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // Constructors, destructors, assignment
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
     StringPhiPsiOutput();
     StringPhiPsiOutput(std::string fname, int n);
     virtual ~StringPhiPsiOutput();


  public:
     virtual void doRun();

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // My Data members
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int numpoints;
};

#endif /*STRINGPHIPSIOUTPUT_H*/
