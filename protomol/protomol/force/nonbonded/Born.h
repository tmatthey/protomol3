#ifndef BORN_H
#define BORN_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/ScalarStructure.h>

namespace ProtoMol {
  class TOneAtomPair;
  class GenericTopology;
  class Vector3DBlock;
  class ScalarStructure;

  template<class TOneAtomPair>
  class Born {
  public:
    void evaluateBorn(TOneAtomPair &myOneAtomPair, GenericTopology *topo,
                      Vector3DBlock *forces, ScalarStructure *energies) {
      // TC: I am commenting in which equations
      //     in my documentation correspond to which terms.
      Real D = 80; // This needs to be not hardcoded once it works
      for (unsigned int i = 0; i < topo->atoms.size(); i++) {
        // IMPORTANT: Eq. 1 assumes redundant pairs for the Sasa Fraction.
        // i.e. (atom1-atom2) and (atom2-atom1) would both contribute
        // ProtoMol does not compute redundant pairs.
        // Therefore, we have to multiply by 2.
        Real f = topo->atoms[i].mySCPISM->sasaFrac;
        Real rw = topo->atoms[i].mySCPISM->R_w;
        Real A_i = topo->atomTypes[topo->atoms[i].type].mySCPISM->A_i;
        Real dR_vdw2 = topo->atoms[i].mySCPISM->dR_vdw2;
        Real gamma = 0.5; // Same for all
        //**********************************************
        // Eq. 3
        // Atoms that are not polar hydrogens have a 0
        // polarFrac, so this will reduce to Eq. 1
        topo->atoms[i].mySCPISM->bornRadius =
          rw + gamma - gamma * dR_vdw2 *
          (A_i + f) + topo->atoms[i].mySCPISM->polarFrac;
        //**********************************************
      }

      for (unsigned int i = 0; i < topo->atoms.size(); i++) {
        // More variable definitions
        Real bornRadius = topo->atoms[i].mySCPISM->bornRadius;
        Real rBornRadius = 1.0 / topo->atoms[i].mySCPISM->bornRadius;
        Real K = (D - 1.0) / 2;
        Real alpha_i = topo->atomTypes[topo->atoms[i].type].mySCPISM->alpha;
        Real Ds_r = (1.0 + D) / (1 + K * exp(-alpha_i * bornRadius)) - 1.0;
        Real rDs_r = 1.0 / Ds_r;
        Real dDs_r = alpha_i * (1.0 / (1.0 + D)) * (1 + Ds_r) * (D - Ds_r);
        Real q_i = topo->atoms[i].scaledCharge;

        for (int k = 0;
             k < myOneAtomPair.getNonbondedForceFunction()->getSize(i); k++) {
          // J is an atom that atom I interacted with at some point
          // during the pairwise computation.
          int j =
            myOneAtomPair.getNonbondedForceFunction()->getInteraction(i, k);

          // I term
          Real forceval = 0.5 * q_i * q_i * rBornRadius * rBornRadius *
            (1.0 - rDs_r - bornRadius * rDs_r * rDs_r * dDs_r) *
            myOneAtomPair.getNonbondedForceFunction()->getDR(i, j);

          // J term
          Real bornRadius_j = topo->atoms[j].mySCPISM->bornRadius;
          Real rBornRadius_j = 1.0 / bornRadius_j;
          Real alpha_j = topo->atomTypes[topo->atoms[j].type].mySCPISM->alpha;
          Real Ds_r_j =
            (1.0 + D) / (1 + K * exp(-alpha_j * bornRadius_j)) - 1.0;
          Real rDs_r_j = 1.0 / Ds_r_j;
          Real dDs_r_j =
            alpha_j * (1.0 / (1.0 + D)) * (1 + Ds_r_j) * (D - Ds_r_j);
          Real q_j = topo->atoms[j].scaledCharge;

          // Add the J term to the force value
          forceval += 0.5 * q_j * q_j * rBornRadius_j * rBornRadius_j *
            (1.0 - rDs_r_j - bornRadius_j * rDs_r_j * rDs_r_j * dDs_r_j) *
            myOneAtomPair.getNonbondedForceFunction()->getDR(j, i);

          //**********************************************
          // Eq. 5
          // Apply the force to atoms i and j.
          // Note: for this we need the diff (x,y,z) and
          // the distance r_ij between them, these are stored
          // in diffs and dists respectively.
          //**********************************************

          (*forces)[i] +=
            myOneAtomPair.getNonbondedForceFunction()->getDiff(i, k) *
            (1.0 / (myOneAtomPair.getNonbondedForceFunction()->getDist(i, k))) *
            forceval;

          (*forces)[j] -=
            myOneAtomPair.getNonbondedForceFunction()->getDiff(i, k) *
            (1.0 / (myOneAtomPair.getNonbondedForceFunction()->getDist(i, k))) *
            forceval;

          //**********************************************
        }

        //**********************************************
        // Eq. 4
        (*energies)[ScalarStructure::COULOMB] += 0.5 * q_i * q_i * rBornRadius *
          (rDs_r - 1);
        //**********************************************
      }

      // Now zero out the fractions.
      for (unsigned int b = 0; b < topo->atoms.size(); b++) {
        topo->atoms[b].mySCPISM->sasaFrac = 0.0;
        topo->atoms[b].mySCPISM->polarFrac = 0.0;
      }

      myOneAtomPair.getNonbondedForceFunction()->clear();
    }
  };
}
#endif
