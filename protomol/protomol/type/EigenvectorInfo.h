/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <protomol/type/Real.h>
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
    EigenvectorInfo() {};
    EigenvectorInfo(unsigned int n, unsigned int m) {
      myEigenvectorLength = n;
      myNumEigenvectors = m;
      myMaxEigenvalue = 0.0;
      myEigenvectors = new Real[n * m * 3];
    }
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
  };
}
#endif /* EIGENVECTORINFO_H */
