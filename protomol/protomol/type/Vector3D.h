/*  -*- c++ -*-  */
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>

namespace ProtoMol {
  //_________________________________________________________________ Vector3D
  /**
   * Container to hold 3D vector/coordinate
   */
  struct Vector3D {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    Real x, y, z;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    Vector3D() : x(0.0), y(0.0), z(0.0) {}

    Vector3D(Real X, Real Y, Real Z) : x(X), y(Y), z(Z) {}

    Vector3D(const Vector3D & c) {
      x = c.x;
      y = c.y;
      z = c.z;
    }

    Vector3D &operator=(const Vector3D &c) {
      x = c.x;
      y = c.y;
      z = c.z;
      return *this;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Vector3D
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:
    /// Index access
    Real operator[](int index) const {
      if (index < 0 || index > 2)
        Report::report << Report::error <<
        "[Vector3D::operator[] const] index out of range" << Report::endr;
      return index == 0 ? x : (index == 1 ? y : z);
    }
    /// Index access
    Real &operator[](int index) {
      if (index < 0 || index > 2)
        Report::report << Report::error <<
        "[Vector3D::operator[] const] index out of range" << Report::endr;
      return index == 0 ? x : (index == 1 ? y : z);
    }

    // Binary operators
    Vector3D operator+(const Vector3D &b) const {
      return Vector3D(x + b.x, y + b.y, z + b.z);
    }
    Vector3D add(const Vector3D &b) const {
      return (*this) + b;
    }

    Vector3D operator-(const Vector3D &b) const {
      return Vector3D(x - b.x, y - b.y, z - b.z);
    }
    Vector3D subtract(const Vector3D &b) const {
      return (*this) - b;
    }

    /// dot product
    Real operator*(const Vector3D &b) const {
      return x * b.x + y * b.y + z * b.z;
    }
    /// dot product
    Real dot(const Vector3D &b) const {
      return (*this) * b;
    }

    /// cross product
    Vector3D operator^(const Vector3D &b) const {
      return Vector3D(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    /// cross product
    Vector3D cross(const Vector3D &b) const {
      return (*this) ^ b;
    }


    Vector3D operator*(Real w) const {
      return Vector3D(x * w, y * w, z * w);
    }
    Vector3D multiply(Real w) const {
      return (*this) * w;
    }

    Vector3D operator/(Real w) const {
      return Vector3D(x / w, y / w, z / w);
    }
    Vector3D divide(Real w) const {
      return (*this) / w;
    }


    // Unary operators
    Vector3D operator-() const {
      return Vector3D(-x, -y, -z);
    }


    // Comparison
    bool operator==(const Vector3D &b) const {
      return x == b.x && y == b.y && z == b.z;
    }

    bool operator!=(const Vector3D &b) const {
      return x != b.x || y != b.y || z != b.z;
    }

    // Assignment
    Vector3D &operator+=(const Vector3D &b) {
      x += b.x;
      y += b.y;
      z += b.z;
      return *this;
    }
    Vector3D &intoAdd(const Vector3D &b) {
      return (*this) += b;
    }

    Vector3D &operator-=(const Vector3D &b) {
      x -= b.x;
      y -= b.y;
      z -= b.z;
      return *this;
    }
    Vector3D &intoSubtract(const Vector3D &b) {
      return (*this) -= b;
    }

    Vector3D &operator*=(Real w) {
      x *= w;
      y *= w;
      z *= w;
      return *this;
    }
    Vector3D &intoMultiply(Real w) {
      return (*this) *= w;
    }

    Vector3D &operator/=(Real w) {
      x /= w;
      y /= w;
      z /= w;
      return *this;
    }
    Vector3D &intoDivide(Real w) {
      return (*this) /= w;
    }


    Vector3D &intoWeightedAdd(Real w, const Vector3D &b) {
      x += w * b.x;
      y += w * b.y;
      z += w * b.z;
      return *this;
    }

    Vector3D &intoWeightedSubtract(Real w, const Vector3D &b) {
      x -= w * b.x;
      y -= w * b.y;
      z -= w * b.z;
      return *this;
    }



    Real normSquared() const {
      return x * x + y * y + z * z;
    }

    Real norm() const {
      return sqrt(x * x + y * y + z * z);
    }

    /// Normalize the Vector3D and return the original length.
    Real normalize() {
      Real len = norm();

      if (len == 0.0)
        return 0.0;

      Real d = 1 / len;
      x *= d; y *= d; z *= d;
      return len;
    }

    /// Return a normalized Vector3D, leave the original Vector3D unchanged.
    Vector3D normalized() const {
      Real len = norm();

      if (len == 0.0) {
        Report::report << Report::recoverable <<
        "[Vector3D::normalized] length is zero." << Report::endr;
        return Vector3D(0.0, 0.0, 0.0);
      }

      return Vector3D(*this / len);
    }

    friend std::ostream &operator<<(std::ostream &OS, const Vector3D &coords) {
      OS << "(" << coords.x << "," << coords.y << "," << coords.z << ")";
      return OS;
    }

    friend std::istream &operator>>(std::istream &OS, Vector3D &coords) {
      OS >> coords.x >> coords.y >> coords.z;
      return OS;
    }

    friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                          const Vector3D &coords) {
      OS << "(" << coords.x << "," << coords.y << "," << coords.z << ")";
      return OS;
    }
  };
}

#endif

