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
    EigenvectorInfo() : myEigenvectorLength ( 0 ), myNumEigenvectors ( 0 ), 
        myEigenvectors ( 0 ), mySingleEigs( 0 ),  
        myEigVecChanged( true ), myMinimumLimit( 0.5 ), currentMode( -1 ) {};

    EigenvectorInfo( unsigned int n, unsigned int m ) : myEigenvectorLength( n ), myNumEigenvectors( m ),
        myMaxEigenvalue( 0.0 ), myEigenvectors ( new double[n * m * 3] ), mySingleEigs( 0 ),
        myEigVecChanged( true ), myMinimumLimit( 0.5 ), currentMode( -1 ) {}

    ~EigenvectorInfo() {
      if ( myEigenvectors ) {
        delete [] myEigenvectors;
        myEigenvectors = 0;
      }
      if ( mySingleEigs ) {
        delete [] mySingleEigs;
        mySingleEigs = 0;
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

    float* getFloatEigPointer() {

      //no eigenvectors assigned?
      if(!myEigenvectors) return (float*)0;

      const unsigned int arrayLen = myEigenvectorLength * myNumEigenvectors * 3;

      //assign storage if required
      if(!mySingleEigs) {  

        try {
          mySingleEigs = new float[arrayLen];
        } catch ( std::bad_alloc& ) {
          return (float*)0;
        }

      }

      //update values if double array updated
      if(myEigVecChanged) {

        for( unsigned int i=0; i<arrayLen; i++)
          mySingleEigs[i] = (float)myEigenvectors[i];

        myEigVecChanged = false;

      }

      return mySingleEigs;

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Eigenvector information
    unsigned int myEigenvectorLength;
    unsigned int myNumEigenvectors;
    double myMaxEigenvalue;
    double *myEigenvectors;

    //OpenMM single precision interface
    float *mySingleEigs;
    bool myEigVecChanged;

    double myMinimumLimit;

    //Analytic integrator
    std::vector< double > myEigenvalues;
    int currentMode;


  };
}
#endif /* EIGENVECTORINFO_H */
