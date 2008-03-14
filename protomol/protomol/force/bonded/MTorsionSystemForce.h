/* -*- c++ -*- */
#ifndef MTORSIONSYSTEMFORCE_H
#define MTORSIONSYSTEMFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/Torsion.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  //____ MTorsionSystemForce
  template<class TBoundaryConditions>
  class MTorsionSystemForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MTorsionSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void calcTorsion(const TBoundaryConditions &boundary,
                     const Torsion &currentTorsion,
                     const Vector3DBlock *positions,
                     Vector3DBlock *forces, Real &energy,
                     ScalarStructure *energies);
    Real calcTorsionEnergy(const TBoundaryConditions &boundary,
                           const Torsion &currentTorsion,
                           const Vector3DBlock *positions);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
  };

  //____ INLINES
  template<class TBoundaryConditions>
  inline void MTorsionSystemForce<TBoundaryConditions>::
  calcTorsion(const TBoundaryConditions &boundary, const Torsion &currTorsion,
              const Vector3DBlock *positions, Vector3DBlock *forces,
              Real &energy, ScalarStructure *energies) {

    int a1 = currTorsion.atom1;
    int a2 = currTorsion.atom2;
    int a3 = currTorsion.atom3;
    int a4 = currTorsion.atom4;

    Vector3D r12(boundary.minimalDifference((*positions)[a2],
                                            (*positions)[a1]));
    Vector3D r23(boundary.minimalDifference((*positions)[a3],
                                            (*positions)[a2]));
    Vector3D r34(boundary.minimalDifference((*positions)[a4],
                                            (*positions)[a3]));

    // Cross product of r12 and r23, represents the plane shared by these two 
    // vectors
    Vector3D a(r12.cross(r23));
    // Cross product of r23 and r34, represents the plane shared by these two 
    // vectors
    Vector3D b(r23.cross(r34));
    // Cross product of r23 and A, represents the plane shared by these two 
    // vectors
    Vector3D c(r23.cross(a));

    // 1/length of Vector A, B and C
    Real ra = 1.0 / a.norm();
    Real rb = 1.0 / b.norm();
    Real rc = 1.0 / c.norm();

    // Normalize A,B and C
    a *= ra;
    b *= rb;
    c *= rc;

    // Calculate phi
    Real cosPhi = a.dot(b);
    Real sinPhi = c.dot(b);
    Real phi = -atan2(sinPhi, cosPhi);

    Real dpotdphi = 0.;
    for (int i = 0; i < currTorsion.multiplicity; i++)

      if (currTorsion.periodicity[i] > 0) {
        dpotdphi -= currTorsion.periodicity[i]
                    * currTorsion.forceConstant[i]
                    * sin(currTorsion.periodicity[i] * phi
                          + currTorsion.phaseShift[i]);

        // Add energy
        energy += currTorsion.forceConstant[i] *
          (1.0 + cos(currTorsion.periodicity[i] * phi +
                     currTorsion.phaseShift[i]));
      } else {
        Real diff = phi - currTorsion.phaseShift[i];

        if (diff < -M_PI)
          diff += 2 * M_PI;
        else if (diff > M_PI)
          diff -= 2 * M_PI;

        dpotdphi += 2.0 * currTorsion.forceConstant[i] * diff;

        // Add energy
        energy += currTorsion.forceConstant[i] * diff * diff;
      }

    // To prevent potential singularities, if abs(sinPhi) <= 0.1, then
    // use another method of calculating the gradient.
    Vector3D f1, f2, f3;
    if (fabs(sinPhi) > 0.1) {
      //  use the sin version to avoid 1/cos terms

      Vector3D dcosdA((a * cosPhi - b) * ra);
      Vector3D dcosdB((b * cosPhi - a) * rb);

      Real k1 = dpotdphi / sinPhi;

      f1.x = k1 * (r23.y * dcosdA.z - r23.z * dcosdA.y);
      f1.y = k1 * (r23.z * dcosdA.x - r23.x * dcosdA.z);
      f1.z = k1 * (r23.x * dcosdA.y - r23.y * dcosdA.x);

      f3.x = k1 * (r23.z * dcosdB.y - r23.y * dcosdB.z);
      f3.y = k1 * (r23.x * dcosdB.z - r23.z * dcosdB.x);
      f3.z = k1 * (r23.y * dcosdB.x - r23.x * dcosdB.y);

      f2.x = k1 *
             (r12.z * dcosdA.y - r12.y * dcosdA.z + r34.y * dcosdB.z - r34.z *
              dcosdB.y);
      f2.y = k1 *
             (r12.x * dcosdA.z - r12.z * dcosdA.x + r34.z * dcosdB.x - r34.x *
              dcosdB.z);
      f2.z = k1 *
             (r12.y * dcosdA.x - r12.x * dcosdA.y + r34.x * dcosdB.y - r34.y *
              dcosdB.x);
    } else {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms

      Vector3D dsindC((c * sinPhi - b) * rc);
      Vector3D dsindB((b * sinPhi - c) * rb);

      Real k1 = -dpotdphi / cosPhi;

      f1.x = k1 *
             ((r23.y * r23.y + r23.z *
               r23.z) * dsindC.x - r23.x * r23.y * dsindC.y - r23.x * r23.z *
              dsindC.z);
      f1.y = k1 *
             ((r23.z * r23.z + r23.x *
               r23.x) * dsindC.y - r23.y * r23.z * dsindC.z - r23.y * r23.x *
              dsindC.x);
      f1.z = k1 *
             ((r23.x * r23.x + r23.y *
               r23.y) * dsindC.z - r23.z * r23.x * dsindC.x - r23.z * r23.y *
              dsindC.y);

      f3 = dsindB.cross(r23) * k1;

      f2.x = k1 *
             (-(r23.y * r12.y + r23.z *
                r12.z) * dsindC.x +
              (2.0 * r23.x * r12.y - r12.x * r23.y) * dsindC.y
              + (2.0 * r23.x * r12.z - r12.x *
                 r23.z) * dsindC.z + dsindB.z * r34.y - dsindB.y * r34.z);
      f2.y = k1 *
             (-(r23.z * r12.z + r23.x *
                r12.x) * dsindC.y +
              (2.0 * r23.y * r12.z - r12.y * r23.z) * dsindC.z
              + (2.0 * r23.y * r12.x - r12.y *
                 r23.x) * dsindC.x + dsindB.x * r34.z - dsindB.z * r34.x);
      f2.z = k1 *
             (-(r23.x * r12.x + r23.y *
                r12.y) * dsindC.z +
              (2.0 * r23.z * r12.x - r12.z * r23.x) * dsindC.x
              + (2.0 * r23.z * r12.y - r12.z *
                 r23.y) * dsindC.y + dsindB.y * r34.x - dsindB.x * r34.y);
    }
    (*forces)[a1] += f1;
    (*forces)[a2] += f2 - f1;
    (*forces)[a3] += f3 - f2;
    (*forces)[a4] -= f3;

    // Add virial
    if (energies->virial()) {
      Real xy = f1.x * r12.y + f2.x * r23.y + f3.x * r34.y;
      Real xz = f1.x * r12.z + f2.x * r23.z + f3.x * r34.z;
      Real yz = f1.y * r12.z + f2.y * r23.z + f3.y * r34.z;

      (*energies)[ScalarStructure::VIRIALXX] += f1.x * r12.x + f2.x * r23.x +
                                                f3.x * r34.x;
      (*energies)[ScalarStructure::VIRIALXY] += xy;
      (*energies)[ScalarStructure::VIRIALXZ] += xz;
      (*energies)[ScalarStructure::VIRIALYX] += xy;
      (*energies)[ScalarStructure::VIRIALYY] += f1.y * r12.y + f2.y * r23.y +
                                                f3.y * r34.y;
      (*energies)[ScalarStructure::VIRIALYZ] += yz;
      (*energies)[ScalarStructure::VIRIALZX] += xz;
      (*energies)[ScalarStructure::VIRIALZY] += yz;
      (*energies)[ScalarStructure::VIRIALZZ] += f1.z * r12.z + f2.z * r23.z +
                                                f3.z * r34.z;
    }
  }

  template<class TBoundaryConditions>
  inline Real MTorsionSystemForce<TBoundaryConditions>::
  calcTorsionEnergy(const TBoundaryConditions &boundary,
                    const Torsion &currTorsion,
                    const Vector3DBlock *positions) {
    int a1 = currTorsion.atom1;
    int a2 = currTorsion.atom2;
    int a3 = currTorsion.atom3;
    int a4 = currTorsion.atom4;

    Vector3D r12 = boundary.minimalDifference((*positions)[a2],
                                              (*positions)[a1]);
    Vector3D r23 = boundary.minimalDifference((*positions)[a3],
                                              (*positions)[a2]);
    Vector3D r34 = boundary.minimalDifference((*positions)[a4],
                                              (*positions)[a3]);

    // Cross product of r12 and r23, represents the plane shared by these two 
    // vectors
    Vector3D a = r12.cross(r23);
    // Cross product of r12 and r23, represents the plane shared by these two
    // vectors
    Vector3D b = r23.cross(r34);

    Vector3D c = r23.cross(a);

    // Calculate phi.
    Real cosPhi = a.dot(b) / (a.norm() * b.norm());
    Real sinPhi = c.dot(b) / (c.norm() * b.norm());
    Real phi = -atan2(sinPhi, cosPhi);

    // Calculate energy.

    Real energy = 0.0;

    for (int i = 0; i < currTorsion.multiplicity; i++)

      if (currTorsion.periodicity[i] > 0)

        // Add energy
        energy += currTorsion.forceConstant[i] *
          (1.0 + cos(currTorsion.periodicity[i] * phi +
                     currTorsion.phaseShift[i]));


      else {
        Real diff = phi - currTorsion.phaseShift[i];

        if (diff < -M_PI) diff += 2 * M_PI;
        else if (diff > M_PI) diff -= 2 * M_PI;

        // Add energy
        energy += currTorsion.forceConstant[i] * diff * diff;
      }

    return energy;
  }
}
#endif /* MTORSIONSYSTEMFORCE_H */
