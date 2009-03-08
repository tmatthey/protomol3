/* -*- c++ -*- */
#ifndef MYOUTPUT_H
#define MYOUTPUT_H

#include <string>
#include <fstream>

class MyOutput {

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // Constructors, destructors, assignment
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
     MyOutput():myfilename("") { myData = NULL; }
     MyOutput(std::string fname): myfilename(fname) { myData = NULL; }
     virtual ~MyOutput() { 
       if (myfile.is_open()) myfile.close();
     }
    
  public:
     virtual void open() {
         if (myfile.is_open()) return;
         else {
            myfile.open(myfilename.c_str(), std::ios::out | std::ios::trunc);
         }
     }

     virtual void open(std::string fname) {
         myfilename = fname;
         if (myfile.is_open()) return;
         else {
            myfile.open(myfilename.c_str(), std::ios::out | std::ios::trunc);
         }
     }
               

  public:
     //
     // The idea : double *d is to be written to the output file.
     // Exactly what data to be written and how will be figured out
     // in derived class.
     //
     virtual void run(double *d){
        myData = d;
        doRun();
     }

     virtual void doRun()=0; //pure virtual function

  protected:
     std::string myfilename;
     std::ofstream myfile;
     double *myData;

};

#endif /*MYOUTPUT_H*/
