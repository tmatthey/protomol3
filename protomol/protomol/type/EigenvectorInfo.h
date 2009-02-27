/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <string>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
namespace ProtoMol
{
  //_________________________________________________________________XYZ
  /**
   * Container holding coordinates/Vector3D and names
   */
  struct EigenvectorInfo {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EigenvectorInfo() : myEigenvectors ( 0 ), currentMode( -1 ), myMinimumLimit( 0.5 ) {};

    EigenvectorInfo( unsigned int n, unsigned int m ) : myEigenvectorLength( n ), myNumEigenvectors( m ),
        myMaxEigenvalue( 0.0 ), myEigenvectors ( new double[n * m * 3] ), currentMode( -1 ), myMinimumLimit( 0.5 ) {}

    ~EigenvectorInfo() {
      if ( myEigenvectors ) {
        delete [] myEigenvectors;
        myEigenvectors = 0;
      }
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class EigenvectorInfo
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool initializeEigenvectors() {
      try {

        myEigenvectors = new double[myEigenvectorLength * myNumEigenvectors * 3];

      } catch ( std::bad_alloc& ) {

        return false;

      }

      return true;

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Eigenvector information
    unsigned int myEigenvectorLength;
    unsigned int myNumEigenvectors;
    double myMaxEigenvalue;
    double *myEigenvectors;

    std::vector< double > myEigenvalues;
    int currentMode;

    double myMinimumLimit;

  };
}
#endif /* EIGENVECTORINFO_H */
