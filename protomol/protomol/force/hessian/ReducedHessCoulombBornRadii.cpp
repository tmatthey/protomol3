#include <protomol/force/hessian/ReducedHessCoulombBornRadii.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/force/coulomb/CoulombBornRadiiForce.h>

using namespace ProtoMol::Report;

using namespace ProtoMol;
//____ ReducedHessCoulombBornRadii
Matrix3By3 ReducedHessCoulombBornRadii::operator()(
  const Real rawEnergy, const Real rawForce, Real a, Real /*rDistSquared*/,
  const Vector3D &rij, const GenericTopology *
  topo, int atom1, int atom2, const Real switchingValue,
  const Real switchingDeriv, const Matrix3By3 &switchingHess,
  ExclusionClass excl,
  const Vector3DBlock *positions, CoulombBornRadiiForce &hForce
  ) const {
  Real D = 80;

  Real na = sqrt(a);

  Real P1r = 0;
  Real P2r = 0;
  //Real tm1 = 1;

  Real f, fPrime, fDoublePrime;
  f = hForce.BornSwitchValue(na);
  fPrime = hForce.BornSwitchDerivative(na);
  fDoublePrime = hForce.BornSwitchSecondDerivative(na);

  if (!topo->atoms[atom1].mySCPISM || !topo->atoms[atom2].mySCPISM ||
      !topo->atomTypes[topo->atoms[atom1].type].mySCPISM ||
      !topo->atomTypes[topo->atoms[atom2].type].mySCPISM)
    report << error << "[ReducedHessCoulombBornRadii::operator()] SCPISM "
      "parameters not set for one or more atoms or atom types." << endr;

  Real bornRadius_i = topo->atoms[atom1].mySCPISM->bornRadius;
  Real bornRadius_j = topo->atoms[atom2].mySCPISM->bornRadius;

  Real K = (D - 1.0) / 2;

  Real sqrt_alpha_i =
    topo->atomTypes[topo->atoms[atom1].type].mySCPISM->sqrt_alpha;
  Real alpha_i = sqrt_alpha_i * sqrt_alpha_i;
  Real Ds_r_i = (1.0 + D) / (1 + K * exp(-alpha_i * bornRadius_i)) - 1.0;
  Real dDs_r_i = alpha_i * (1.0 / (1.0 + D)) * (1 + Ds_r_i) * (D - Ds_r_i);
  //second derivative
  Real d2Ds_r2_i = (alpha_i / (1 + D)) * (D - 1 - 2 * Ds_r_i) * dDs_r_i;

  Real gamma = 0.5;
  int type_i = topo->atoms[atom1].type;
  Real B_i = topo->atomTypes[type_i].mySCPISM->B_i;
  Real C_i = topo->atomTypes[type_i].mySCPISM->C_i;
  Real dR_vdw2_i = topo->atoms[atom1].mySCPISM->dR_vdw2;

  Real dR_i_coeff = gamma * B_i * dR_vdw2_i * (1 / na) *
    (fPrime * exp(-C_i * na) - C_i * f * exp(-C_i * na));

  Real sqrt_alpha_j =
    topo->atomTypes[topo->atoms[atom2].type].mySCPISM->sqrt_alpha;
  Real alpha_j = sqrt_alpha_j * sqrt_alpha_j;
  Real Ds_r_j = (1.0 + D) / (1 + K * exp(-alpha_j * bornRadius_j)) - 1.0;
  Real dDs_r_j = alpha_j * (1.0 / (1.0 + D)) * (1 + Ds_r_j) * (D - Ds_r_j);
  //second derivative
  Real d2Ds_r2_j = (alpha_j / (1 + D)) * (D - 1 - 2 * Ds_r_j) * dDs_r_j;

  int type_j = topo->atoms[atom2].type;
  Real B_j = topo->atomTypes[type_j].mySCPISM->B_i;
  Real C_j = topo->atomTypes[type_j].mySCPISM->C_i;
  Real dR_vdw2_j = topo->atoms[atom2].mySCPISM->dR_vdw2;

  Real dR_j_coeff = gamma * B_j * dR_vdw2_j * (1 / na) * 
    (fPrime * exp(-C_j * na) - C_j * f * exp(-C_j * na));

  Real expCi = exp(-C_i * na);
  Real expCj = exp(-C_j * na);

  Real d2R_dr_i_r_coeff =
    gamma * B_i * dR_vdw2_i *
    ((1 / (na * na)) * (fDoublePrime * expCi - 2 * C_i * fPrime * expCi +
                        C_i * C_i * f * expCi) - (1 / (na * na * na)) *
     (fPrime * expCi - C_i * f * expCi));

  Real d2R_dr_j_r_coeff =
    gamma * B_j * dR_vdw2_j *
    ((1 / (na * na)) * (fDoublePrime * expCj - 2 * C_j * fPrime * expCj + C_j *
                        C_j * f * expCj) - (1 / (na * na * na)) *
     (fPrime * expCj - C_j * f * expCj));

  //
  //Find out if there is polar interaction between atom1 and atom2
  //

  //
  //Between atom1 and atom2
  //

  if (topo->atomTypes[topo->atoms[atom1].type].mySCPISM->isHbonded == PH &&
      topo->atomTypes[topo->atoms[atom2].type].mySCPISM->isHbonded == PA &&
      (topo->atoms[atom1].residue_seq !=
       topo->atoms[atom2].residue_seq)) {   // Polar
    Real E_i = 0.80;     // Currently all polar H+ have this value
    Real g_i;
    if (topo->atoms[atom1].name == "HN" && topo->atoms[atom2].name == "O")
      g_i = -0.378;
    else g_i = topo->atomTypes[type_i].mySCPISM->g_i;
    Real g_j = topo->atomTypes[type_j].mySCPISM->g_i;

    dR_i_coeff += g_i * g_j * (1 / na) * (fPrime * exp(-E_i * na) + f * 
                                          exp(-E_i * na) * -E_i);

    d2R_dr_i_r_coeff +=
      g_i * g_j * ((1 / (na * na)) * (fDoublePrime * exp(-E_i * na) - 2 * E_i *
                                      fPrime * exp(-E_i * na) + E_i * E_i * f *
                                      exp(-E_i * na)) - (1 / (na * na * na)) *
                   (fPrime * exp(-E_i * na) + f * exp(-E_i * na) * -E_i));
  }

  //
  //Between atom2 and atom1
  //
  if (topo->atomTypes[type_j].mySCPISM->isHbonded == PH &&
      topo->atomTypes[type_i].mySCPISM->isHbonded == PA &&
      (topo->atoms[atom2].residue_seq !=
       topo->atoms[atom1].residue_seq)) {
    Real E_j = 0.80;         // Currently all polar H+ have this value
    Real g_i;
    if (topo->atoms[atom2].name == "HN" &&
        topo->atoms[atom1].name == "O")
      g_i = -0.378;
    else
      g_i = topo->atomTypes[type_j].mySCPISM->g_i;
    Real g_j = topo->atomTypes[type_i].mySCPISM->g_i;

    dR_j_coeff += g_i * g_j * (1 / na) *
      (fPrime * exp(-E_j * na) + f * exp(-E_j * na) * -E_j);

    d2R_dr_j_r_coeff += g_i * g_j *
      ((1 / (na * na)) * (fDoublePrime * exp(-E_j * na) - 2 * E_j * fPrime *
                          exp(-E_j * na) + E_j * E_j * f * exp(-E_j * na)) -
       (1 / (na * na * na)) * 
       (fPrime * exp(-E_j * na) + f * exp(-E_j * na) * -E_j));
  }

  Real R_i = topo->atoms[atom1].mySCPISM->bornRadius;
  Real R_j = topo->atoms[atom2].mySCPISM->bornRadius;

  Real MR_i = (1 / (R_i * R_i)) - (1 / (R_i * R_i * Ds_r_i)) -
    (1 / (R_i * Ds_r_i * Ds_r_i)) * dDs_r_i;
  Real MR_j = (1 / (R_j * R_j)) - 
    (1 / (R_j * R_j * Ds_r_j)) - (1 / (R_j * Ds_r_j * Ds_r_j)) * dDs_r_j;

  Real dMR_i_dR = (2 / (Ds_r_i * R_i * R_i * R_i)) -
    (2 / (R_i * R_i * R_i)) + (2 / (R_i * R_i * Ds_r_i * Ds_r_i)) * dDs_r_i +
    (2 / (R_i * Ds_r_i * Ds_r_i * Ds_r_i)) * dDs_r_i * dDs_r_i -
    (1 / (Ds_r_i * Ds_r_i * R_i)) * d2Ds_r2_i;

  Real dMR_j_dR = (2 / (Ds_r_j * R_j * R_j * R_j)) -
    (2 / (R_j * R_j * R_j)) + (2 / (R_j * R_j * Ds_r_j * Ds_r_j)) * dDs_r_j +
    (2 / (R_j * Ds_r_j * Ds_r_j * Ds_r_j)) * dDs_r_j * dDs_r_j -
    (1 / (Ds_r_j * Ds_r_j * R_j)) * d2Ds_r2_j;

  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_j = topo->atoms[atom2].scaledCharge;

  P1r = 0.5 * scaledCharge_i * scaledCharge_i * MR_i * dR_i_coeff + 0.5 *
        scaledCharge_j * scaledCharge_j * MR_j * dR_j_coeff;

  P2r = 0.5 * scaledCharge_i * scaledCharge_i *
        (dMR_i_dR * dR_i_coeff * dR_i_coeff + MR_i *
         d2R_dr_i_r_coeff) + 0.5 * scaledCharge_j * scaledCharge_j *
        (dMR_j_dR * dR_j_coeff * dR_j_coeff + MR_j * d2R_dr_j_r_coeff);

  Matrix3By3 vec_rij_ij(rij, rij);   // outer products of vectors r_ij
  Matrix3By3 I(1, 0, 0, 0, 1, 0, 0, 0, 1);   // now I is identity matrix
  Matrix3By3 H((I * P1r + vec_rij_ij * P2r));

  for (unsigned int ii = 0; ii < topo->atoms.size(); ii++) {
    if (((int)ii == atom1) || ((int)ii == atom2)) continue;
    Vector3D rk =
      topo->minimalDifference((*positions)[atom1], (*positions)[ii]);

    Real rik = rk.norm();
    f = hForce.BornSwitchValue(rik);
    fPrime = hForce.BornSwitchDerivative(rik);
    fDoublePrime = hForce.BornSwitchSecondDerivative(rik);

    Real dR_isum_coeff = gamma * B_i * dR_vdw2_i * (1 / rik) *
      (fPrime * exp(-C_i * rik) - C_i * f * exp(-C_i * rik));

    Matrix3By3 tk(rk, rij);
    Real fac_i = 0.5 * scaledCharge_i * scaledCharge_i *
      (dMR_i_dR * dR_isum_coeff * dR_i_coeff);
    H += tk * fac_i;
  }

  for (unsigned int ii = 0; ii < topo->atoms.size(); ii++) {
    if (((int)ii == atom1) || ((int)ii == atom2)) continue;
    Vector3D rk =
      topo->minimalDifference((*positions)[atom2], (*positions)[ii]);

    Real rjk = rk.norm();
    f = hForce.BornSwitchValue(rjk);
    fPrime = hForce.BornSwitchDerivative(rjk);
    fDoublePrime = hForce.BornSwitchSecondDerivative(rjk);

    Real dR_jsum_coeff = gamma * B_j * dR_vdw2_j * (1 / rjk) *
      (fPrime * exp(-C_j * rjk) - C_j * f * exp(-C_j * rjk));

    Matrix3By3 tjk(rk, -rij);
    Real fac_j = 0.5 * scaledCharge_i * scaledCharge_i *
      (dMR_j_dR * dR_jsum_coeff * dR_j_coeff);
    H += tjk * fac_j;
  }

  for (unsigned int ii = 0; ii < topo->atoms.size(); ii++) {
    if (((int)ii == atom1) || ((int)ii == atom2)) continue;

    Real scaledCharge_ii = topo->atoms[ii].scaledCharge;
    //cout<<"Atom "<<ii<<", charge "<<scaledCharge_ii<<endl;
    int type_ii = topo->atoms[ii].type;
    Real B_ii = topo->atomTypes[type_ii].mySCPISM->B_i;
    Real sqrt_alpha_ii =
      topo->atomTypes[topo->atoms[ii].type].mySCPISM->sqrt_alpha;
    Real alpha_ii = sqrt_alpha_ii * sqrt_alpha_ii;
    Real dr_vdw2_ii = topo->atoms[ii].mySCPISM->dR_vdw2;
    Real C_ii = topo->atomTypes[type_ii].mySCPISM->C_i;
    Real R_ii = topo->atoms[ii].mySCPISM->bornRadius;
    Real Ds_r_ii = (1 + D) / (1 + K * exp(-alpha_ii * R_ii)) - 1;
    Real dDs_r_ii =
      alpha_ii * (1.0 / (1.0 + D)) * (1 + Ds_r_ii) * (D - Ds_r_ii);
    Real d2Ds_r2_ii = (alpha_ii / (1 + D)) * (D - 1 - 2 * Ds_r_ii) * dDs_r_ii;

    Vector3D rn1 =
      topo->minimalDifference((*positions)[ii], (*positions)[atom1]);
    Vector3D rn2 =
      topo->minimalDifference((*positions)[ii], (*positions)[atom2]);
    Real ria1 = rn1.norm();
    //Real ria12=ria1*ria1;
    f = hForce.BornSwitchValue(ria1);
    fPrime = hForce.BornSwitchDerivative(ria1);
    fDoublePrime = hForce.BornSwitchSecondDerivative(ria1);

    Real dR_ii_coeff_1 = gamma * B_ii * dr_vdw2_ii * (1 / ria1) * 
      (fPrime * exp(-C_ii * ria1) - C_ii * f * exp(-C_ii * ria1));
    if (topo->atomTypes[topo->atoms[ii].type].mySCPISM->isHbonded == PH &&
        topo->atomTypes[topo->atoms[atom1].type].mySCPISM->isHbonded == PA &&
        (topo->atoms[atom1].residue_seq !=
         topo->atoms[ii].residue_seq)) {   // Polar
      Real E_i = 0.80;   // Currently all polar H+ have this value
      Real g_i;
      if (topo->atoms[ii].name == "HN" && topo->atoms[atom1].name == "O")
        g_i = -0.378;
      else g_i = topo->atomTypes[topo->atoms[ii].type].mySCPISM->g_i;
      Real g_j = topo->atomTypes[topo->atoms[atom1].type].mySCPISM->g_i;

      dR_ii_coeff_1 += g_i * g_j * (1 / ria1) *
        (fPrime * exp(-E_i * ria1) - f * E_i * exp(-E_i * ria1));
    }

    Real dMR_ii_dR =
      (2 / (Ds_r_ii * R_ii * R_ii * R_ii)) - (2 / (R_ii * R_ii * R_ii)) +
      (2 / (R_ii * R_ii * Ds_r_ii * Ds_r_ii)) * dDs_r_ii +
      (2 / (R_ii * Ds_r_ii * Ds_r_ii * Ds_r_ii)) * dDs_r_ii * dDs_r_ii -
      (1 / (Ds_r_ii * Ds_r_ii * R_ii)) * d2Ds_r2_ii;

    Real ria2 = rn2.norm();
    //Real ria22=ria2*ria2;
    f = hForce.BornSwitchValue(ria2);
    fPrime = hForce.BornSwitchDerivative(ria2);
    fDoublePrime = hForce.BornSwitchSecondDerivative(ria2);
    Real dR_ii_coeff_2 = gamma * B_ii * dr_vdw2_ii * (1 / ria2) *
      (fPrime * exp(-C_ii * ria2) - C_ii * f * exp(-C_ii * ria2));

    if (topo->atomTypes[topo->atoms[ii].type].mySCPISM->isHbonded == PH &&
        topo->atomTypes[topo->atoms[atom2].type].mySCPISM->isHbonded == PA &&
        (topo->atoms[atom2].residue_seq !=
         topo->atoms[ii].residue_seq)) {   // Polar
      Real E_i = 0.80;   // Currently all polar H+ have this value
      Real g_i;
      if (topo->atoms[ii].name == "HN" && topo->atoms[atom2].name == "O")
        g_i = -0.378;
      else g_i = topo->atomTypes[topo->atoms[ii].type].mySCPISM->g_i;
      Real g_j = topo->atomTypes[topo->atoms[atom2].type].mySCPISM->g_i;

      dR_ii_coeff_2 += g_i * g_j * (1 / ria2) *
        (fPrime * exp(-E_i * ria2) - f * E_i * exp(-E_i * ria2));
    }

    Matrix3By3 vec_rin1_rin2(rn1, rn2);
    Matrix3By3 temp = vec_rin1_rin2 * dMR_ii_dR * dR_ii_coeff_1 *
      dR_ii_coeff_2 * 0.5 * scaledCharge_ii * scaledCharge_ii;
    H -= temp;
  }

  return H;
}

