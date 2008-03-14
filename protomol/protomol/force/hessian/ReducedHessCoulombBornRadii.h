/*  -*- c++ -*-  */
#ifndef REDUCEDHESSCOULOMBBORNRADII_H
#define REDUCEDHESSCOULOMBBORNRADII_H

#include <protomol/type/Matrix3By3.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  class GenericTopology;
  class CoulombBornRadiiForce;

  class ReducedHessCoulombBornRadii {
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  Constructor,destructor
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ReducedHessCoulombBornRadii
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Matrix3By3 operator()(const Real rawEnergy, const Real rawForce,
                          Real distSquared, Real rDistSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2, const Real switchingValue,
                          const Real switchingDeriv,
                          const Matrix3By3 &switchingHess,
                          ExclusionClass excl, const Vector3DBlock *positions,
                          CoulombBornRadiiForce &hForce) const;
  };
}
#endif
