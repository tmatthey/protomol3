/*  -*- c++ -*-  */
#ifndef ATOM_H
#define ATOM_H

#include <protomol/type/Real.h>
#include <string>
#include <vector>

namespace ProtoMol {
  //________________________________________ Atom

  struct SCPISMAtomParameters {
    SCPISMAtomParameters() : sasaFrac(0.0), polarFrac(0.0) {}
    Real sqrtalphaSCPISM;
    Real alphaSCPISM;
    ///< For SCPISM calculations
    ///< Parameterized by Atom Type according to Hassan et al. (2002)
    Real sasaFrac;
    ///< For SCPISM calculations

    Real bornRadius, D_s;
    ///< For SCPISM calculations

    Real polarFrac;
    ///< For SCPISM calculations

    Real dR_vdw2;
    Real r_cov;
    Real R_iw;
    Real R_w;
    Real R_p;
    //##CRS
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
             deltaM(0.0), deltaQ(0.0), Qold(0.0), Qnew(0.0), stageNumber(0),
             mySCPISM(0) {}

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

    Real deltaM, deltaQ;
    ///< For iSG simulations.  This difference between the mass and charge,
    ///< respectively, of the atom's two identities.

    Real Qold, Qnew;
    ///< For iSG simulations.  The atom's charge for its old and new identities.
    ///< Needed to correctly compute DeltaMu for intramolecular interactions.

    int stageNumber;
    ///< For iSG simulations.  Molecule transformations can be broken into 
    ///< stages so that only small fragments of the whole molecule can be 
    ///< transformed at any time.  This is critical to being able to transform 
    ///< a small molecule into a relatively large molecule.  This variable is 
    ///< the number of the stage in which this atom will be transformed.  For 
    ///< example, suppose I break a molecular transformation into 5 stages.  
    ///< If this atom's stageNumber = 3, that means that this atom will have
    ///< its identity transformed only when 2 <= lambda <= 3.  -- TIM

    Real alphaLJ;
    ///< For iSG simulations.  Atoms which are dummy atoms in one particular 
    ///< identity but are real, interacting atoms in another identity require 
    ///< a positive, non-zero alpha parameter for the soft-core Lennard-Jones 
    ///< function used by iSGMD.  Atoms which are real, interacting atoms in 
    ///< both identities of a transformation attempt need to have alphaLJ = 0 
    ///< for the soft-core LJ function.  I found that when transforming a water
    ///< oxygen atom into an alcohol oxygen atom using the standard alphaLJ
    ///< value of 0.5, the simulation was wildly unstable, but if I used 
    ///< alphaLJ = 0 only for the oxygen atom then the transformation proceeded
    ///< smoothly.  -- TIM 3/31/2005  (alphaLJ will be set prior to force
    ///< calculations by ModifierISG::pickNewMolecule(...) )

    std::string name;
    std::string residue_name;
    int residue_seq;
    SCPISMAtomParameters *mySCPISM;
  };
}
#endif /* ATOM_H */

