/*  -*- c++ -*-  */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <protomol/type/Real.h>
#include <string>

namespace ProtoMol {
  /**
   * Repository namespace of constants.
   */
  namespace Constant {
    extern const std::string PROTOMOL_HR;
    extern const std::string PRINTINDENT;
    extern const unsigned int PRINTMAXWIDTH;

    extern const int FASTDELTAMAX;
    extern const Real SQRTCOULOMBCONSTANT;
    extern const Real PRESSUREFACTOR;
    extern const Real BOLTZMANN;
    extern const Real PDBVELSCALINGFACTOR;

    extern const Real MAXREAL;
    extern const Real MINREAL;
    extern const Real REAL_INFINITY;
    extern const Real REAL_NAN;
    extern const int MAX_INT;
    extern const int MAX_INT_2;
    extern const Real EPSILON;
    extern const Real TINY;

    extern const Real TIMEFACTOR;
    extern const Real INV_TIMEFACTOR;
    extern const Real PERIODIC_BOUNDARY_TOLERANCE;

    extern const Real EPS_GOURAUD_THRESHOLD;
    extern const Real EPS_SMOOTH_LINE_FACTOR;

    namespace SI {
      extern const Real C;
      extern const Real COULOMB_FACTOR;
      extern const Real ELECTRON_CHARGE;
      extern const Real LENGTH_AA;
      extern const Real AVOGADRO;
      extern const Real AMU;
      extern const Real KCAL;
      extern const Real TIME_FS;
      extern const Real BOLTZMANN;
    }
    /// Returns actual limits/values of the numerical constants
    std::string print();
  }
}

#endif
