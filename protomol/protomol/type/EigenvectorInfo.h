/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <protomol/type/Real.h>
#include <protomol/type/Vector3DBlock.h>
#include <string>
#include <iostream>

using std::cout;
using std::endl;
namespace ProtoMol {
  //_________________________________________________________________XYZ
  /**
   * Container holding coordinates/Vector3D and names
   */
  struct EigenvectorInfo {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EigenvectorInfo() {
      myEigenvectors = 0;
      currentMode = -1;
    };

    EigenvectorInfo(unsigned int n, unsigned int m) {
      myEigenvectorLength = n;
      myNumEigenvectors = m;
      myMaxEigenvalue = 0.0;
      myEigenvectors = new Real[n * m * 3];      
      currentMode = -1;
    }

    ~EigenvectorInfo() {
      if(myEigenvectors){
        delete [] myEigenvectors;
        myEigenvectors = 0;
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class EigenvectorInfo
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void initializeEigenvectors() {
      //myEigenvectors = new Real[myEigenvectorLength*myNumEigenvectors];
      myEigenvectors = new Real[myEigenvectorLength * myNumEigenvectors * 3];
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Eigenvector information
    unsigned int myEigenvectorLength;
    unsigned int myNumEigenvectors;
    Real myMaxEigenvalue;
    Real *myEigenvectors;
    Vector3DBlock myEigenvalues;
    int currentMode;
  };
}
#endif /* EIGENVECTORINFO_H */
