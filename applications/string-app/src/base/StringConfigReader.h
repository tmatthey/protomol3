/* -*- c++ -*- */
#ifndef STRINGCONFIGREADER_H
#define STRINGCONFIGREADER_H

#include <protomol/io/Reader.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace ProtoMol;

class StringConfigReader :
     public Reader {

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // Constructors, destructors, assignment
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   public:
       StringConfigReader();
       explicit StringConfigReader (const std::string &filename);
       virtual ~StringConfigReader();
       

       std::vector<std::string> *GetConfigFilenames() { return &configfiles; }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   public:
       virtual bool tryFormat();
       virtual bool read();


   
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //  data members 
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   public :
     std::vector<std::string> configfiles;
      

};

#endif /*STRINGCONFIGREADER_H*/
