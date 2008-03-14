/* -*- c++ -*- */
#ifndef PERODICBOUNDARYCONDITIONS_H
#define PERODICBOUNDARYCONDITIONS_H

#include <string>

#include <protomol/config/Parameter.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  //________________________________________ PeriodicBoundaryConditions

  /**
   * Implements periodic boundary conditions, defining how we measure distances
   * and accounting the wrapping-around effect.
   * The class use a couple of shorts cut's to avoid rint and to many div's and
   * mul's.
   */
  class PeriodicBoundaryConditions {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    PeriodicBoundaryConditions();
    PeriodicBoundaryConditions(const Vector3D &e1, const Vector3D &e2,
                               const Vector3D &e3, const Vector3D &origin);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PeriodicBoundaryConditions
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Set method for the dimensions of the (original) simulation box.
    void set(const Vector3D &e1, const Vector3D &e2, const Vector3D &e3,
             const Vector3D &origin);

    /// Perform a minimal-image subtraction.
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2) const {
      Vector3D diff(c2);
      diff -= c1;
      if (myOrthogonal)
        if (diff.normSquared() > myD) {
          // diff not small engouh, we do have to do more ...
          // ... may be it's just a single wrapping
          if (diff.x < -myH2.x)
            diff.x += myH.x;
          else if (diff.x > myH2.x)
            diff.x -= myH.x;
          if (diff.y < -myH2.y)
            diff.y += myH.y;
          else if (diff.y > myH2.y)
            diff.y -= myH.y;
          if (diff.z < -myH2.z)
            diff.z += myH.z;
          else if (diff.z > myH2.z)
            diff.z -= myH.z;
          if (diff.normSquared() > myD) {
            // ... distance was pretty big, hu? ...
            diff.x -= myE1.x * rint(myE1r.x * diff.x);
            diff.y -= myE2.y * rint(myE2r.y * diff.y);
            diff.z -= myE3.z * rint(myE3r.z * diff.z);
          }
        }
      else
        // ... really difficult!
        diff -=
          Vector3D(myE1 * rint(myE1r.dot(diff)) + myE2 * rint(myE2r.dot(
                diff)) + myE3 * rint(myE3r.dot(diff)));
      return diff;
    }

    // Perform a minimal-image subtraction and computes the squared distance
    Vector3D minimalDifference(const Vector3D &c1, const Vector3D &c2,
                               Real &distSquared) const {
      Vector3D diff(c2);
      diff -= c1;
      if (myOrthogonal) {
        distSquared = diff.normSquared();
        if (distSquared > myD) {
          // diff not small engouh, we do have to do more ...
          // ... may be it's just a single wrapping
          if (diff.x < -myH2.x)
            diff.x += myH.x;
          else if (diff.x > myH2.x)
            diff.x -= myH.x;
          if (diff.y < -myH2.y)
            diff.y += myH.y;
          else if (diff.y > myH2.y)
            diff.y -= myH.y;
          if (diff.z < -myH2.z)
            diff.z += myH.z;
          else if (diff.z > myH2.z)
            diff.z -= myH.z;
          distSquared = diff.normSquared();
          if (distSquared > myD) {
            // ... distance was pretty big, hu? ...
            diff.x -= myE1.x * rint(myE1r.x * diff.x);
            diff.y -= myE2.y * rint(myE2r.y * diff.y);
            diff.z -= myE3.z * rint(myE3r.z * diff.z);
            distSquared = diff.normSquared();
          }
        }
      } else {
        // ... really difficult!
        diff -=
          Vector3D(myE1 * rint(myE1r.dot(diff)) + myE2 * rint(myE2r.dot(
                diff)) + myE3 * rint(myE3r.dot(diff)));
        distSquared = diff.normSquared();
      }
      return diff;
    }

    /// Find the position in the basis/original cell/image.
    Vector3D minimalPosition(const Vector3D &c) const {
      Vector3D diff(c);
      diff -= myOrigin;
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return diff;
        else {
          // ... may be it's just a single wrapping
          if (diff.x < -myH2.x)
            diff.x += myH.x;
          else if (diff.x > myH2.x)
            diff.x -= myH.x;
          if (diff.y < -myH2.y)
            diff.y += myH.y;
          else if (diff.y > myH2.y)
            diff.y -= myH.y;
          if (diff.z < -myH2.z)
            diff.z += myH.z;
          else if (diff.z > myH2.z)
            diff.z -= myH.z;
          if (diff.normSquared() <= myD)
            return diff;
          else
            // ... distance was pretty big, hu? ...
            return Vector3D(diff - Vector3D(myE1.x * rint(myE1r.x * diff.x),
                myE2.y * rint(myE2r.y * diff.y),
                myE3.z * rint(myE3r.z * diff.z)));
        }
      else
        // ... really difficult!
        return Vector3D(diff - myE1 * rint(myE1r.dot(diff))
          - myE2 * rint(myE2r.dot(diff))
          - myE3 * rint(myE3r.dot(diff)));
    }

    /// Find the lattice vector difference between two positions
    Vector3D minimalTranslationDifference(const Vector3D &c1,
                                          const Vector3D &c2) const {
      Vector3D diff(c2 - c1);
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return Vector3D(0.0, 0.0, 0.0);
        else
          return Vector3D(myE1.x * rint(myE1r.x * diff.x),
            myE2.y * rint(myE2r.y * diff.y),
            myE3.z * rint(myE3r.z * diff.z));
      else
        // ... really difficult!
        return Vector3D(myE1 * rint(myE1r.dot(diff))
          + myE2 * rint(myE2r.dot(diff))
          + myE3 * rint(myE3r.dot(diff)));
    }

    /// Find the lattice translation relative to the original cell/image
    Vector3D minimalTranslationPosition(const Vector3D &c) const {
      Vector3D diff(c - myOrigin);
      if (myOrthogonal)
        if (diff.normSquared() <= myD)
          // diff so small, we do not have to do more ...
          return Vector3D(0.0, 0.0, 0.0);
        else
          return Vector3D(myE1.x * rint(myE1r.x * diff.x),
            myE2.y * rint(myE2r.y * diff.y),
            myE3.z * rint(myE3r.z * diff.z));
      else
        // ... really difficult!
        return Vector3D(myE1 * rint(myE1r.dot(diff))
          + myE2 * rint(myE2r.dot(diff))
          + myE3 * rint(myE3r.dot(diff)));

    }

    /// basis/unit vector e1
    const Vector3D &e1()     const {return myE1;}
    /// basis/unit vector e2
    const Vector3D &e2()     const {return myE2;}
    /// basis/unit vector e3
    const Vector3D &e3()     const {return myE3;}
    /// inverse basis/unit vector e1
    const Vector3D &e1r()    const {return myE1r;}
    /// inverse basis/unit vector e2
    const Vector3D &e2r()    const {return myE2r;}
    /// inverse basis/unit vector e3
    const Vector3D &e3r()    const {return myE3r;}
    /// origin of the minimal image
    const Vector3D &origin() const {return myOrigin;}
    /// minimal corner of the bounding box of the minimal image/cell
    const Vector3D &getMin() const {return myMin;}
    /// maximal corner of the bounding box of the minimal image/cell
    const Vector3D &getMax() const {return myMax;}
    Real getVolume()         const {return myV;};
    bool isOrthogonal()      const {return myOrthogonal;};

    /// Boolean's defining boundary conditions
    enum {PERIODIC = 1, VACUUM = 0};

    /// Builds the lattice vector for a given cutoff. (0,0,0) is not included.
    std::vector<Vector3D> buildLatticeVectors(Real cutoff) const;

    /// Returns the keyword of the boundary conditions
    const std::string &getKeyword() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameters) const;
    static PeriodicBoundaryConditions make(std::vector<Value> values);
    static unsigned int getParameterSize() {return 4;}

    /// Returns possible default values for the parameters based on the
    /// positions
    std::vector<Parameter> getDefaults(const Vector3DBlock &positions) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3D myE1;
    Vector3D myE2;
    Vector3D myE3;
    Vector3D myE1r;
    Vector3D myE2r;
    Vector3D myE3r;
    Vector3D myOrigin;
    Vector3D myMin;
    Vector3D myMax;

    Real myDX;
    Real myDY;
    Real myDZ;
    /// maximal distance between two positions where plain subtraction if safe
    Real myD;  
    Vector3D myH;
    Vector3D myH2;

    Real myV;
    bool myOrthogonal;
  };
}
#endif /* PERODICBOUNDARYCONDITIONS_H */
