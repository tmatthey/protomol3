/*  -*- c++ -*-  */
#ifndef ATOM_H
#define ATOM_H

#include <protomol/type/Real.h>
#include <string>
#include <vector>

namespace ProtoMol {
  //________________________________________ Atom

  struct SCPISMAtomParameters {
    SCPISMAtomParameters() {}

    //Pre force initialization
    void preForce(){
      bornRadius = zeta;
      energySum = true;
      D_s = 0.0;
    }

    Real sqrtalphaSCPISM;

    ///< For SCPISM calculations
    //Updated version CRS 01/24/09

    Real bornRadius, D_s;
    Real zeta, eta;
    bool energySum;
  };
  /**
   * This class defines the information for one atom.  It contains the type
   * of the atom, its charge (scaled by a constant factor), and the next
   * atom in this atom's cell list.  All other pieces of atom information
   * are either in the atom type or the coordinate vectors (for positions,
   * velocities, and so on).
   */
  struct Atom {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Atom() : type(-1), scaledCharge(0.0), scaledMass(0.0), hvyAtom(-1),
      atomNum(-1), cellListNext(-1), molecule(-1), mybonds(std::vector<int>()),
             mySCPISM_A(0) {}

    // A default constructor for the atom class

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int type;
    ///< The type number of the atom (an index into the array of atomTypes).

    Real scaledCharge;
    ///< The charge of this atom, scaled by the square root of the Coulomb
    ///< constant.

    Real scaledMass;
    ///< The mass of this atom (can be modified, e.g., iSG)

    int hvyAtom;
    ///< The size of the Group for Heavy Atom Come First ordering
    ///< for hydrogen hvyAtom=0, else, hvyAtom = 1+#of attached hydrogens

    int atomNum;
    ///< Original order number of the atom, used to undo Heavy Atom Come
    ///< First ordering

    mutable int cellListNext;
    ///< The index of the next atom in this atom's cell list, or -1 if this
    ///< atom is the last in its list.

    int molecule;
    ///< The ID# of the molecule to which this atom belongs

    std::vector<int> mybonds;
    ///< this vector holds the bond IDs for those bonds made with this atom

    std::string name;
    std::string residue_name;
    int residue_seq;
    SCPISMAtomParameters *mySCPISM_A;
  };
}
#endif /* ATOM_H */

