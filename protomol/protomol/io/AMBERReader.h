/*  -*- c++ -*-  */
#ifndef AMBERREADER_H
#define AMBERREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/PSF.h>
#include <protomol/type/PAR.h>

namespace ProtoMol {

  //_________________________________________________________________PSFReader
  class AMBERReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    AMBERReader(); 
    explicit AMBERReader(const std::string& filename);
    virtual ~AMBERReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open(){return File::open();}
    virtual bool open(const std::string& filename){return File::open(filename);}
    virtual bool open(const char* filename){return File::open(filename);}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PSF
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(PSF& psf, PAR& par);
    //bool read(PSF psf, PAR par);
    //bool readAllAmber(PSF& psf, PAR& par);

    PSF* orphanPSF();
    PAR* orphanPAR();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend AMBERReader& operator>>(AMBERReader& topReader, PSF& psf);
    friend AMBERReader& operator>>(AMBERReader& topReader, PAR& par);
    
    void writeData();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //private:
    //AMBER* myAMBER;
  public:
    PSF *myPSF;
    PAR *myPAR;
    PAR::CharmmTypeEnum myCharmmType;
    PAR::CharmmTypeEnum myCharmmTypeDetected;

  };

  //____________________________________________________________________________INLINES

}
#endif /* AMBERREADER_H */
