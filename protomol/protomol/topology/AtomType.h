/*  -*- c++ -*-  */
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <string>

#include <protomol/type/Real.h>

namespace ProtoMol {
  enum Hbonded {NO = 0, PH = 1, PA = 2};

  struct SCPISMAtomTypeParameters {
    Real sqrt_alpha;        ///< JI: sqrt of alpha i used in CoulombSCPISM
    Real alpha;             ///< JI: alpha i used in CoulombSCPISM
    Real g_i;               ///< TC: used for polar H+, in addendum
    Hbonded isHbonded;      ///< TC: type of H-bond, none, polar, or acceptor
    Real A_i;               ///< JI: A_i in (2) on SCPISM.impl
    Real B_i;               ///< JI: B_i in (2) on SCPISM.impl
    Real C_i;               ///< JI: C_i in (2) on SCPISM.impl
  };

  //________________________________________ AtomType
  /**
   * This class contains information common to one type of atom.  This
   * currently only includes the name and mass of the atom type.  The
   * Lennard-Jones parameters are stored in the LennardJonesParameterTable
   * structure.
   */
  struct AtomType {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AtomType() : name(""), mass(0.0), charge(0.0), symbolName(""), mySCPISM(0)
    {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::string name;       ///< The name of this atom type.
    Real mass;              ///< The mass of this atom type.
    Real charge;            ///< The charge of this atom type.
    std::string symbolName; ///< The symbol untity name of this atom type.
    SCPISMAtomTypeParameters *mySCPISM;
  };
}
#endif /* ATOMTYPE_H */
