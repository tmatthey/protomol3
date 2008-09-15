#include <protomol/integrator/hessian/Hessian.h>

#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/force/hessian/ReducedHessAngle.h>
#include <protomol/force/LennardJonesForce.h>
#include <protomol/force/CoulombForce.h>
#include <protomol/switch/CnSwitchingFunction.h>
#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/switch/C1SwitchingFunction.h>

#include <protomol/force/hessian/HessDihedral.h>
#include <protomol/force/hessian/ReducedHessBond.h>
#include <protomol/force/hessian/ReducedHessCoulomb.h>
#include <protomol/force/hessian/ReducedHessCoulombDiElec.h>
#include <protomol/force/hessian/ReducedHessCoulombSCPISM.h>
#include <protomol/force/hessian/ReducedHessCoulombBornRadii.h>
#include <protomol/force/hessian/ReducedHessLennardJones.h>
#include <protomol/force/coulomb/CoulombForceDiElec.h>
#include <protomol/force/coulomb/CoulombSCPISMForce.h>
#include <protomol/force/coulomb/CoulombBornRadiiForce.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//constructors
Hessian::Hessian() {
  hessM = 0;
}

Hessian::Hessian(unsigned int szin) {
  sz = szin;
  try{
    hessM = new double[sz * sz];   //assign array
  }catch(bad_alloc&){
    report << error << "[Hessian::Hessian] Cannot allocate memory for Hessian!" << endr;
  }
}

Hessian::~Hessian() {
  if (hessM != 0) delete[] hessM;
}

// copy constructor
Hessian::Hessian(const Hessian &hess) {
  sz = hess.sz;
  if (hess.hessM != 0) {
    try{
        hessM = new double[sz * sz];   //assign array
    }catch(bad_alloc&){
        report << error << "[Hessian::Hessian] Cannot allocate memory for Hessian!" << endr;
    }
    //hessM = new double[sz * sz];
    for (unsigned int i = 0; i < sz * sz; i++) hessM[i] = hess.hessM[i];
  } else
    hessM = 0;
  myBond = hess.myBond;
  myAngle = hess.myAngle;
  myCoulomb = hess.myCoulomb;
  myCoulombDielec = hess.myCoulombDielec;
  myCoulombSCPISM = hess.myCoulombSCPISM;
  myCoulombBornRadii = hess.myCoulombBornRadii;
  myLennardJones = hess.myLennardJones;
  myDihedral = hess.myDihedral;
  myImproper = hess.myImproper;
  cCutoff = hess.cCutoff;
  cSwitchon = hess.cSwitchon;
  cSwitch = hess.cSwitch;
  cOrder = hess.cOrder;
  cSwitchoff = hess.cSwitchoff;
  lCutoff = hess.lCutoff;
  lSwitchon = hess.lSwitchon;
  lSwitch = hess.lSwitch;
  lOrder = hess.lOrder;
  lSwitchoff = hess.lSwitchoff;
}

void Hessian::initialData(unsigned int szin) {
  sz = szin;
  if (hessM == 0){
      try{
        hessM = new double[sz * sz];   //assign array
      }catch(bad_alloc&){
        report << error << "[Hessian::initialData] Cannot allocate memory for Hessian!" << endr;
      }
  }
  //if (hessM == 0) hessM = new double[sz * sz]; //assign array
}

void Hessian::findForces(ForceGroup *overloadedForces) {
  vector<Force *> ListForces = overloadedForces->getForces();
  //
  lCutoff = lSwitchon = lSwitch = cCutoff = cSwitchon = cSwitch = 0.0;
  lOrder = cOrder = lSwitchoff = cSwitchoff = 0.0;
  D = 78.0; S = 0.3; epsi = 1.0;
  myBond = myAngle = myCoulomb = myCoulombDielec = myCoulombSCPISM =
                                                     myCoulombBornRadii =
                                                       myLennardJones =
                                                         myDihedral =
                                                           myImproper = false;
  for (unsigned int i = 0; i < ListForces.size(); i++) {
    if (equalNocase(ListForces[i]->getId(), "Bond"))
      myBond = true;
    else if (equalNocase(ListForces[i]->getId(), "Angle"))
      myAngle = true;
    else if (equalStartNocase("CoulombDiElec", ListForces[i]->getId())) {
      myCoulombDielec = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text))
            cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-D"))
          D = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-S"))
          S = Fparam[j].value;
        else if (equalNocase(Fparam[j].keyword, "-EPS"))
          epsi = Fparam[j].value;
      }
    } else if (equalStartNocase("CoulombSCPISM", ListForces[i]->getId())) {
      myCoulombSCPISM = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("CoulombBornRadii", ListForces[i]->getId())) {
      myCoulombBornRadii = true;
      swt = 1;         //default
    } else if (equalStartNocase("Coulomb", ListForces[i]->getId())) {
      myCoulomb = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          cCutoff = Fparam[j].value;
          if (cSwitch == 0) cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          cSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) cSwitch = 3;
          else cSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          cOrder = Fparam[j].value;
          cSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          cSwitchoff = Fparam[j].value;
          cSwitch = 3;
        }
      }
    } else if (equalStartNocase("LennardJones", ListForces[i]->getId())) {
      myLennardJones = true;
      vector<Parameter> Fparam;
      ListForces[i]->getParameters(Fparam);
      for (unsigned int j = 0; j < Fparam.size(); j++) {
        if (equalNocase(Fparam[j].keyword, "-cutoff")) {
          lCutoff = Fparam[j].value;
          if (lSwitch == 0) lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-switchon")) {
          lSwitchon = Fparam[j].value;
          if (equalStartNocase("C2", Fparam[j].text)) lSwitch = 2;
          else if (equalStartNocase("Cn", Fparam[j].text)) lSwitch = 3;
          else lSwitch = 1;
        } else if (equalNocase(Fparam[j].keyword, "-n")) {
          lOrder = Fparam[j].value;
          lSwitch = 3;
        } else if (equalNocase(Fparam[j].keyword, "-switchoff")) {
          lSwitchoff = Fparam[j].value;
          lSwitch = 3;
        }
      }
    } else if (equalStartNocase("Dihedral", ListForces[i]->getId()))
      myDihedral = true;
    else if (equalStartNocase("Improper", ListForces[i]->getId()))
      myImproper = true;
  }
}

void Hessian::evaluate(const Vector3DBlock *myPositions,
                       const GenericTopology *myTopo,
                       bool mrw) {
  int a1, a2, a3, eye, jay;
  Real tempf, ms1, ms2, ms3;
  unsigned int i;
  ReducedHessAngle rh;
  Matrix3By3 rha;
  ExclusionClass ec;

  sz = 3 * myPositions->size();
  //Impropers
  if (myImproper) {
    HessDihedral hi;        //create improper hessian
    for (unsigned int i = 0; i < myTopo->impropers.size(); i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->impropers[i].multiplicity; j++)
        if (myTopo->impropers[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = myTopo->impropers[i].atom1; aout[1] =
          myTopo->impropers[i].atom2;
        aout[2] = myTopo->impropers[i].atom3; aout[3] =
          myTopo->impropers[i].atom4;
        //
        hi.evaluate(myTopo->impropers[i], (*myPositions)[aout[0]],
                    (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                    (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++)
          for (int kk = 0; kk < 4; kk++) {
            int mi = ii * 36 + kk * 3;            //kk * 36 + ii*3;
            Matrix3By3 rhd(hi.hessD[mi], hi.hessD[mi + 1], hi.hessD[mi + 2],
                           hi.hessD[mi + 12], hi.hessD[mi + 13],
                           hi.hessD[mi + 14],
                           hi.hessD[mi + 24], hi.hessD[mi + 25],
                           hi.hessD[mi + 26]);
            for (int ll = 0; ll < 3; ll++)
              for (int mm = 0; mm < 3; mm++)
                if (mrw)
                  hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                    ((rhd(ll, mm) / sqrt(myTopo->atoms[aout[ii]].scaledMass *
                                         myTopo->atoms[aout[kk]].scaledMass)));
                else hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                    rhd(ll, mm);

          }

      }
    }
  }
  //Dihedrals
  if (myDihedral) {
    HessDihedral hd;        //create dihedral hessian
    for (unsigned int i = 0; i < myTopo->dihedrals.size(); i++) {
      bool nonZForce = false;       //test for force constants
      for (int j = 0; j < myTopo->dihedrals[i].multiplicity; j++)
        if (myTopo->dihedrals[i].forceConstant[j]) nonZForce = true;

      if (nonZForce) {
        int aout[4];
        aout[0] = myTopo->dihedrals[i].atom1; aout[1] =
          myTopo->dihedrals[i].atom2;
        aout[2] = myTopo->dihedrals[i].atom3; aout[3] =
          myTopo->dihedrals[i].atom4;
        //
        hd.evaluate(myTopo->dihedrals[i], (*myPositions)[aout[0]],
                    (*myPositions)[aout[1]], (*myPositions)[aout[2]],
                    (*myPositions)[aout[3]]);
        //output sparse matrix
        for (int ii = 0; ii < 4; ii++)
          for (int kk = 0; kk < 4; kk++) {
            int mi = ii * 36 + kk * 3;            //kk * 36 + ii*3;
            Matrix3By3 rhd(hd.hessD[mi], hd.hessD[mi + 1], hd.hessD[mi + 2],
                           hd.hessD[mi + 12], hd.hessD[mi + 13],
                           hd.hessD[mi + 14],
                           hd.hessD[mi + 24], hd.hessD[mi + 25],
                           hd.hessD[mi + 26]);
            for (int ll = 0; ll < 3; ll++)
              for (int mm = 0; mm < 3; mm++)
                if (mrw)
                  hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                    ((rhd(ll, mm) / sqrt(myTopo->atoms[aout[ii]].scaledMass *
                                         myTopo->atoms[aout[kk]].scaledMass)));
                else hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                    rhd(ll, mm);

          }

      }
    }
  }
  //Bonds
  if (myBond)
    for (i = 0; i < myTopo->bonds.size(); i++) {
      a1 = myTopo->bonds[i].atom1; a2 = myTopo->bonds[i].atom2;
      Real r_0 = myTopo->bonds[i].restLength;
      Real k = myTopo->bonds[i].springConstant;
      Matrix3By3 bondHess12 =
        reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
      //output sparse matrix
      for (int ll = 0; ll < 3; ll++)
        for (int mm = 0; mm < 3; mm++) {
          eye = a1 * 3 + 1; jay = a2 * 3 + 1;
          tempf = bondHess12(ll, mm);
          if (mrw) {
            ms1 = tempf / sqrt(myTopo->atoms[a1].scaledMass *
                               myTopo->atoms[a1].scaledMass);
            ms2 = tempf / sqrt(myTopo->atoms[a2].scaledMass *
                               myTopo->atoms[a2].scaledMass);
            ms3 = -tempf / sqrt(myTopo->atoms[a1].scaledMass *
                                myTopo->atoms[a2].scaledMass);
          } else {
            ms2 = ms1 = tempf;
            ms3 = -tempf;
          }
          hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
          hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
          hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
          hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
        }

    }

  //Angles
  if (myAngle)
    for (i = 0; i < myTopo->angles.size(); i++) {
      a1 = myTopo->angles[i].atom1;
      a2 = myTopo->angles[i].atom2;
      a3 = myTopo->angles[i].atom3;
      Real theta0 = myTopo->angles[i].restAngle;
      Real k_t = myTopo->angles[i].forceConstant;
      Real ubConst = myTopo->angles[i].ureyBradleyConstant;
      Real ubRestL = myTopo->angles[i].ureyBradleyRestLength;
      // ReducedHessAngle for atoms a1, a2 and a3
      rh.evaluate((*myPositions)[a1], (*myPositions)[a2], (*myPositions)[a3],
                  k_t,
                  theta0);
      //ureyBradley
      if (ubConst) {
        //Cheat using bond hessian as same as UB!!!!
        Matrix3By3 ubm =
          reducedHessBond((*myPositions)[a1], (*myPositions)[a3], ubConst,
                          ubRestL);
        rh.accumulateTo(0, 0, ubm);
        rh.accumulateTo(2, 2, ubm);
        rh.accumulateNegTo(2, 0, ubm);
        rh.accumulateNegTo(0, 2, ubm);
      }
      //output sparse matrix
      int aout[3];
      aout[0] = a1; aout[1] = a2; aout[2] = a3;
      for (int ii = 0; ii < 3; ii++)
        for (int kk = 0; kk < 3; kk++) {
          rha = rh(ii, kk);
          for (int ll = 0; ll < 3; ll++)
            for (int mm = 0; mm < 3; mm++)
              if (mrw)
                hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                  ((rha(ll, mm) / sqrt(myTopo->atoms[aout[ii]].scaledMass *
                                       myTopo->atoms[aout[kk]].scaledMass)));
              else hessM[(aout[ii] * 3 + ll) * sz + aout[kk] * 3 + mm] +=
                  rha(ll, mm);

        }

    }

  //Lennard jones
  if (myLennardJones)
    for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
      for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++) {
        //if not bonded/dihedral
        ec = myTopo->exclusions.check(i, j);
        if (ec != EXCLUSION_FULL) {
          Vector3D rij =
            myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
          //Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
          Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          Real a = rij.normSquared();
          ReducedHessLennardJones rHess;
          LennardJonesForce hForce;
          Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
          if (lSwitch) {
            if (lSwitch == 3) {
              CnSwitchingFunction cnsf(lSwitchon, lCutoff, lOrder, lSwitchoff);
              cnsf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = cnsf.hessian(rij, a);
              }
            } else if (lSwitch == 2) {
              C2SwitchingFunction c2sf(lSwitchon, lCutoff);
              c2sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c2sf.hessian(rij, a);
              }
            } else {
              C1SwitchingFunction c1sf(lCutoff);
              c1sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c1sf.hessian(rij, a);
              }
            }
          }
          rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz,
                      ec);
          //output sparse matrix
          for (int ll = 0; ll < 3; ll++)
            for (int mm = 0; mm < 3; mm++) {
              eye = i * 3 + 1; jay = j * 3 + 1;
              tempf = rha(ll, mm);
              if (mrw) {
                ms1 = tempf / sqrt(myTopo->atoms[i].scaledMass *
                                   myTopo->atoms[i].scaledMass);
                ms2 = tempf / sqrt(myTopo->atoms[j].scaledMass *
                                   myTopo->atoms[j].scaledMass);
                ms3 = -tempf / sqrt(myTopo->atoms[i].scaledMass *
                                    myTopo->atoms[j].scaledMass);
              } else {
                ms2 = ms1 = tempf;
                ms3 = -tempf;
              }
              hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
              hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
              hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
              hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
            }

        } else {}
      }

  //Coulombic
  if (myCoulomb)
    for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
      for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++) {
        //if not bonded/dihedral
        ec = myTopo->exclusions.check(i, j);
        if (ec != EXCLUSION_FULL) {
          Vector3D rij =
            myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
          //Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
          Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          Real a = rij.normSquared();
          ReducedHessCoulomb rHess;
          CoulombForce hForce;
          Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
          if (cSwitch) {
            if (cSwitch == 3) {
              CnSwitchingFunction cnsf(cSwitchon, cCutoff, cOrder, cSwitchoff);
              cnsf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = cnsf.hessian(rij, a);
              }
            } else if (cSwitch == 2) {
              C2SwitchingFunction c2sf(cSwitchon, cCutoff);
              c2sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c2sf.hessian(rij, a);
              }
            } else {
              C1SwitchingFunction c1sf(cCutoff);
              c1sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c1sf.hessian(rij, a);
              }
            }
          }
          rha =
            rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec);
          //output sparse matrix
          for (int ll = 0; ll < 3; ll++)
            for (int mm = 0; mm < 3; mm++) {
              eye = i * 3 + 1; jay = j * 3 + 1;
              tempf = rha(ll, mm);
              if (mrw) {
                ms1 = tempf / sqrt(myTopo->atoms[i].scaledMass *
                                   myTopo->atoms[i].scaledMass);
                ms2 = tempf / sqrt(myTopo->atoms[j].scaledMass *
                                   myTopo->atoms[j].scaledMass);
                ms3 = -tempf / sqrt(myTopo->atoms[i].scaledMass *
                                    myTopo->atoms[j].scaledMass);
              } else {
                ms2 = ms1 = tempf;
                ms3 = -tempf;
              }
              hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
              hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
              hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
              hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
            }

        } else {}
      }

  //Coulombic Implicit solvent
  if (myCoulombDielec)
    for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
      for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++) {
        //if not bonded/dihedral
        ec = myTopo->exclusions.check(i, j);
        if (ec != EXCLUSION_FULL) {
          Vector3D rij =
            myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
          //Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
          Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          Real a = rij.normSquared();
          ReducedHessCoulombDielec rHess;
          CoulombForceDiElec hForce;
          Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
          if (cSwitch) {
            if (cSwitch == 3) {
              CnSwitchingFunction cnsf(cSwitchon, cCutoff, cOrder, cSwitchoff);
              cnsf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = cnsf.hessian(rij, a);
              }
            } else if (cSwitch == 2) {
              C2SwitchingFunction c2sf(cSwitchon, cCutoff);
              c2sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c2sf.hessian(rij, a);
              }
            } else {
              C1SwitchingFunction c1sf(cCutoff);
              c1sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                mz = c1sf.hessian(rij, a);
              }
            }
          }
          rha =
            rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec,
                  D, S,
                  epsi);
          //output sparse matrix
          for (int ll = 0; ll < 3; ll++)
            for (int mm = 0; mm < 3; mm++) {
              eye = i * 3 + 1; jay = j * 3 + 1;
              tempf = rha(ll, mm);
              if (mrw) {
                ms1 = tempf / sqrt(myTopo->atoms[i].scaledMass *
                                   myTopo->atoms[i].scaledMass);
                ms2 = tempf / sqrt(myTopo->atoms[j].scaledMass *
                                   myTopo->atoms[j].scaledMass);
                ms3 = -tempf / sqrt(myTopo->atoms[i].scaledMass *
                                    myTopo->atoms[j].scaledMass);
              } else {
                ms2 = ms1 = tempf;
                ms3 = -tempf;
              }
              hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
              hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
              hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
              hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
            }

        } else {}
      }

  if (myCoulombSCPISM)
    evaluateCoulombSCPISM(myPositions, myTopo, mrw);

  if (myCoulombBornRadii)
    evaluateCoulombBornRadii(myPositions, myTopo, mrw);
  //PrintHessian();
}

void Hessian::evaluateCoulombSCPISM(const Vector3DBlock *myPositions,
                                    const GenericTopology *myTopo,
                                    bool mrw) {
  int eye, jay;
  Real tempf, ms1, ms2, ms3;
  Matrix3By3 rha;
  ExclusionClass ec;

  report << plain << "Evaluating Coulomb SCPISM " << endr;

  for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
    for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++) {
      //if not bonded/dihedral
      ec = myTopo->exclusions.check(i, j);
      if (ec != EXCLUSION_FULL) {
        Vector3D rij =
          myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
        //Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
        Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Real a = rij.normSquared();
        ReducedHessCoulombSCPISM rHess;
        CoulombSCPISMForce hForce;
        Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
        if (!myTopo->atoms[i].mySCPISM || !myTopo->atoms[j].mySCPISM)
          report << error <<
          "[Hessian::evaluateCoulombSCPISM] SCPISM data not set." << endr;
        Real alpha_ij = myTopo->atoms[i].mySCPISM->sqrtalphaSCPISM *
                        myTopo->atoms[j].mySCPISM->sqrtalphaSCPISM;
        if (cSwitch) {
          if (cSwitch == 3) {
            CnSwitchingFunction cnsf(cSwitchon, cCutoff, cOrder, cSwitchoff);
            cnsf(swtchV, swtchD, a);
            if (swtchV > 0.0 && swtchV < 1.0) {
              hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
              mz = cnsf.hessian(rij, a);
            }
          } else if (cSwitch == 2) {
            C2SwitchingFunction c2sf(cSwitchon, cCutoff);
            c2sf(swtchV, swtchD, a);
            if (swtchV > 0.0 && swtchV < 1.0) {
              hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
              mz = c2sf.hessian(rij, a);
            }
          } else {
            C1SwitchingFunction c1sf(cCutoff);
            c1sf(swtchV, swtchD, a);
            if (swtchV > 0.0 && swtchV < 1.0) {
              hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
              mz = c1sf.hessian(rij, a);
            }
          }
        }
        rha =
          rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec,
                alpha_ij,
                80.0);
        for (int ll = 0; ll < 3; ll++)
          for (int mm = 0; mm < 3; mm++) {
            eye = i * 3 + 1; jay = j * 3 + 1;
            tempf = rha(ll, mm);
            if (mrw) {
              ms1 = tempf / sqrt(myTopo->atoms[i].scaledMass *
                                 myTopo->atoms[i].scaledMass);
              ms2 = tempf / sqrt(myTopo->atoms[j].scaledMass *
                                 myTopo->atoms[j].scaledMass);
              ms3 = -tempf / sqrt(myTopo->atoms[i].scaledMass *
                                  myTopo->atoms[j].scaledMass);
            } else {
              ms2 = ms1 = tempf;
              ms3 = -tempf;
            }
            hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
            hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
            hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
            hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
          }

      } else {}
    }

}

void Hessian::evaluateCoulombBornRadii(const Vector3DBlock *myPositions,
                                       const GenericTopology *myTopo,
                                       bool mrw) {
  int eye, jay;
  Real tempf, ms1, ms2, ms3;
  Matrix3By3 rha, rhc;
  ExclusionClass ec;

  for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
    for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++) {
      //if not bonded/dihedral
      ec = myTopo->exclusions.check(i, j);
      if (ec != EXCLUSION_FULL) {
        Vector3D rij =
          myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
        //Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
        Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Real a = rij.normSquared();

        ReducedHessCoulombBornRadii rHess;
        swt = myTopo->doSCPISM;
        CoulombBornRadiiForce hForce(swt);
        Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
        rha =
          rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec,
                myPositions,
                hForce);
        for (int ll = 0; ll < 3; ll++)
          for (int mm = 0; mm < 3; mm++) {
            eye = i * 3 + 1; jay = j * 3 + 1;
            tempf = rha(ll, mm);
            if (mrw) {
              ms1 = tempf / sqrt(myTopo->atoms[i].scaledMass *
                                 myTopo->atoms[i].scaledMass);
              ms2 = tempf / sqrt(myTopo->atoms[j].scaledMass *
                                 myTopo->atoms[j].scaledMass);
              ms3 = -tempf / sqrt(myTopo->atoms[i].scaledMass *
                                  myTopo->atoms[j].scaledMass);
            } else {
              ms2 = ms1 = tempf;
              ms3 = -tempf;
            }
            hessM[(eye + ll - 1) * sz + eye + mm - 1] += ms1;
            hessM[(jay + ll - 1) * sz + jay + mm - 1] += ms2;
            hessM[(eye + ll - 1) * sz + jay + mm - 1] += ms3;
            hessM[(jay + ll - 1) * sz + eye + mm - 1] += ms3;
          }

      } else {}
    }

}

void Hessian::clear() {
  if (hessM != 0)
    for (unsigned int i = 0; i < sz * sz; i++) hessM[i] = 0.0;

}

