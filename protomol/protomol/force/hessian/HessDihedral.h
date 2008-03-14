/* -*- c++ -*- */
#ifndef HESSDIHEDRAL_H
#define HESSDIHEDRAL_H

#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  /**
   *
   *  Calculates Hessian for Diherdal PE term.
   *  Rotates the dihedral into a bi-planar configuration
   *  to reduct the number of calculations. Then rotates
   *  the results back again
   */
  class HessDihedral {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    HessDihedral();

    //torsion data including positions etc
    HessDihedral(const Torsion &currTorsion,
                 const Vector3D &a1, const Vector3D &a2,
                 const Vector3D &a3, const Vector3D &a4);
    ~HessDihedral() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HessDihedral
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void evaluate(const Torsion &currTorsion, const Vector3D &a1,
                  const Vector3D &a2, const Vector3D &a3,
                  const Vector3D &a4);

  private:
    //Use aRot to rotate the vector back into real space
    double *rotateV3D(double *aRot, double *mf);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class HessDihedral
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    double hessD[144];  //storage for 4x4 hessian matrices in 3D

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };
}
#endif
