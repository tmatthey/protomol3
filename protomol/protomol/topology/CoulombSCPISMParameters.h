/* -*- c++ -*- */
#ifndef COULOMBSCPISM_H
#define COULOMBSCPISM_H

#include <protomol/type/Real.h>
#include <protomol/topology/AtomType.h>

namespace ProtoMol {
  //________________________________________ CoulombSCPISMParameters

  /// The Coulomb SCPISM for an atom type
  struct CoulombSCPISMParameters {
  public:

    CoulombSCPISMParameters(Real al = 0.0, Real hb = 0.5, Real r = 0.5,
                            Real sq = 0.0, Real rc = 0.37, Real ga = 0.52,
                            Hbonded is = NO)
      : alpha_i(al), hbond_factor(hb), R_iw(r), sqrt_alpha_i(sq),
        r_cov(rc), gamma_i(ga), isHbonded(is) {}
    
    void set(Real al = 0.0, Real hb = 0.5, Real r = 0.5, Real sq = 0.0,
             Real rc = 0.37, Real ga = 0.0052, Hbonded is = NO,
             Real A = 0.0, Real B = 0.0, Real C = 0.0, Real Rvdw = 0.0) {
      alpha_i = al;
      hbond_factor = hb;
      R_iw = r;
      sqrt_alpha_i = sq;
      r_cov = rc;
      gamma_i = ga;
      isHbonded = is;
      A_i = A;
      B_i = B;
      C_i = C;
      R_vdw = Rvdw;
    }

    Real alpha_i; // Alpha_i controls slope of D(r) around atom type i
    Real hbond_factor; // hbond_factor (Polar H) * hbond_factor (PA) controls
                       // increment of Born radius to correct h bonding strength
    Real R_iw; // extension of Born radius R_iw to obtain R_ip
    Real sqrt_alpha_i; // Square root of alpha_i
    Real r_cov; // Covalent radius (only values for C,N,O,S,H)
    Real gamma_i;  // hydrophobic energy term
    Hbonded isHbonded; // whether atom is involved in H-bonding
    Real A_i; // A_i in (2) on SCPISM.DOC
    Real B_i; // B_i in (2) on SCPISM.DOC
    Real C_i; // C_i in (2) on SCPISM.DOC
    Real R_vdw; // R_i,vdw in (2) on SCPISM.DOC
  };
}
#endif /* not COULOMBSCPISM_H */
