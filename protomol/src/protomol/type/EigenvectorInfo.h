/*  -*- c++ -*-  */
#ifndef EIGENVECTORINFO_H
#define EIGENVECTORINFO_H

#include <vector>
#include <protomol/type/Real.h>
namespace ProtoMol {
	/**
	 * Container holding coordinates/Vector3D and names
	 */
	struct EigenvectorInfo {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		EigenvectorInfo();
		EigenvectorInfo( unsigned int n, unsigned int m );

		~EigenvectorInfo();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class EigenvectorInfo
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		bool initializeEigenvectors();
		float *getFloatEigPointer();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Eigenvector information
		unsigned int myEigenvectorLength;
		unsigned int myNumEigenvectors;
		unsigned int myNumUsedEigenvectors;
		Real myMaxEigenvalue;
		Real *myEigenvectors;

		//Current and original max subspace eigenvalue,
		//for adaptive timestep
		Real myOrigCEigval, myNewCEigval;
		Real myOrigTimestep;

		//re-diagonalization flag
		bool reDiagonalize;
		bool havePositionsChanged;
		
		bool OpenMMMinimize;
		unsigned int RediagonalizationCount;

		//OpenMM single precision interface
		float *mySingleEigs;
		bool myEigVecChanged;

		Real myMinimumLimit;

		//Analytic integrator
		std::vector< Real > myEigenvalues;
		int currentMode;
	};
}

#endif // EIGENVECTORINFO_H
