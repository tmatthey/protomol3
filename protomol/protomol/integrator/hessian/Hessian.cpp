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

namespace ProtoMol {
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
  //set maxiumum cutoff value
  cutOff = max(cSwitchoff, lSwitchoff);
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
            Matrix3By3 rhd = hi(ii, kk);
            outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
                myTopo->atoms[aout[kk]].scaledMass, rhd, mrw, sz, hessM);
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
            Matrix3By3 rhd = hd(ii, kk);
            outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
                myTopo->atoms[aout[kk]].scaledMass, rhd, mrw, sz, hessM);
          }

      }
    }
  }

  //Bonds
  if (myBond){
    for (i = 0; i < myTopo->bonds.size(); i++) {
      a1 = myTopo->bonds[i].atom1; a2 = myTopo->bonds[i].atom2;
      Real r_0 = myTopo->bonds[i].restLength;
      Real k = myTopo->bonds[i].springConstant;
      Matrix3By3 bondHess12 =
        reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
      //output sparse matrix
      outputSparsePairMatrix(a1,a2,myTopo->atoms[a1].scaledMass,myTopo->atoms[a2].scaledMass,
                                bondHess12,mrw,sz,hessM);

    }
  }

  //Angles
    if (myAngle){
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
          outputSparseMatrix(aout[ii], aout[kk], myTopo->atoms[aout[ii]].scaledMass,
              myTopo->atoms[aout[kk]].scaledMass, rha, mrw, sz, hessM);

        }

    }
  }
  //Pairwise forces
  for (unsigned int i = 0; i < myTopo->atoms.size(); i++){
    for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++){ 
      Matrix3By3 rhp(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      //Lennard jones
      if (myLennardJones) 
        rhp += evaluatePairsMatrix(i, j, LENNARDJONES, myPositions, myTopo, mrw);
      //Coulombic
      if (myCoulomb)
        rhp += evaluatePairsMatrix(i, j, COULOMB, myPositions, myTopo, mrw);
      //Coulombic Implicit solvent
      if (myCoulombDielec)
        rhp += evaluatePairsMatrix(i, j, COULOMBDIELEC, myPositions, myTopo, mrw);
      //SCP
      if (myCoulombSCPISM)
        rhp += evaluatePairsMatrix(i, j, COULOMBSCPISM, myPositions, myTopo, mrw);
      //Bourn radii
      if (myCoulombBornRadii)  
        evaluateCoulombBornRadiiPair(i, j, myPositions, myTopo, mrw, i, j, sz, hessM);
      //output sum to matrix
      outputSparsePairMatrix(i, j, myTopo->atoms[i].scaledMass, myTopo->atoms[j].scaledMass,
                              rhp, mrw, sz, hessM);
      //
    }
  }
  //
}

void Hessian::evaluatePairs(int i, int j, int pairType, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo, bool mrw, 
                                            int mat_i, int mat_j, int mat_sz, double * mat_array) {

  Matrix3By3 rha = evaluatePairsMatrix(i, j, pairType, myPositions, myTopo, mrw);
  //
  outputSparsePairMatrix(mat_i,mat_j,myTopo->atoms[i].scaledMass,myTopo->atoms[j].scaledMass,
                          rha,mrw,mat_sz,mat_array);

}

Matrix3By3 Hessian::evaluatePairsMatrix(int i, int j, int pairType, const Vector3DBlock *myPositions,
                                        const GenericTopology *myTopo, bool mrw) {

  Matrix3By3 rha(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  ExclusionClass ec;
  //Coulomb
  ReducedHessCoulomb rHess;
  CoulombForce hForce;
  ReducedHessCoulombDielec rHessDi;
  CoulombForceDiElec hForceDi;
  ReducedHessCoulombSCPISM rHessScp;
  CoulombSCPISMForce hForceScp;
  //LJ
  ReducedHessLennardJones rHessLj;
  LennardJonesForce hForceLj;
  //SCP
  Real alpha_ij;
  //Switches
  Real pCutoff, pSwitchon, pSwitch, pOrder, pSwitchoff; 


  if(pairType == LENNARDJONES){	//setup switch paramiters
    pCutoff = lCutoff; pSwitchon = lSwitchon; 
    pSwitch = lSwitch; pOrder = lOrder; pSwitchoff = lSwitchoff; 
  }else{
    pCutoff = cCutoff; pSwitchon = cSwitchon; 
    pSwitch = cSwitch; pOrder = cOrder; pSwitchoff = cSwitchoff; 
  }
  //if not bonded/dihedral
  ec = myTopo->exclusions.check(i, j);
  if (ec != EXCLUSION_FULL) {
      Vector3D rij =
            myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
      //
      Matrix3By3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      Real a = rij.normSquared();
      Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
      //SCP extras
      if(pairType == COULOMBSCPISM){
          if (!myTopo->atoms[i].mySCPISM || !myTopo->atoms[j].mySCPISM)
                    report << error << "[Hessian::evaluateCoulombSCPISM] SCPISM data not set." << endr;
          alpha_ij = myTopo->atoms[i].mySCPISM->sqrtalphaSCPISM *
                            myTopo->atoms[j].mySCPISM->sqrtalphaSCPISM;
      }
      //
      if (pSwitch) {
          if (pSwitch == 3) {
              CnSwitchingFunction cnsf(pSwitchon, pCutoff, pOrder, pSwitchoff);
              cnsf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                  }
                  mz = cnsf.hessian(rij, a);
              }
          } else if (pSwitch == 2) {
              C2SwitchingFunction c2sf(pSwitchon, pCutoff);
              c2sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                break;
                  }
                  mz = c2sf.hessian(rij, a);
              }
          } else {
              C1SwitchingFunction c1sf(pCutoff);
              c1sf(swtchV, swtchD, a);
              if (swtchV > 0.0 && swtchV < 1.0) {
                  switch(pairType){
                    case LENNARDJONES:	hForceLj(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);;
                                        break;
                    case COULOMB:		hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBDIELEC:	hForceDi(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                        break;
                    case COULOMBSCPISM: hForceScp(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
                                break;
                  }
                  mz = c1sf.hessian(rij, a);
              }
          }
      }
      switch(pairType){
          case LENNARDJONES:  rha = rHessLj(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec);
                              break;
          case COULOMB:		rha = rHess(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec);
                              break;
          case COULOMBDIELEC:	rha = rHessDi(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec, D, S, epsi);
                              break;
          case COULOMBSCPISM: rha = rHessScp(rawE, rawF, a, a, rij, myTopo, (unsigned)i, (unsigned)j, swtchV, swtchD, mz, ec, alpha_ij, 80.0);
                              break;
      }
      //
  }
  return rha;
  //
}

void Hessian::outputSparsePairMatrix(int i, int j, Real massi, Real massj, 
                                        Matrix3By3 rha, bool mrw, int arrSz, double *basePoint){
    int eye, jay;
    Real tempf, ms1, ms2, ms3;

    //output sparse matrix
    for (int ll = 0; ll < 3; ll++){
        for (int mm = 0; mm < 3; mm++) {
            eye = i * 3 + 1; jay = j * 3 + 1;
            tempf = rha(ll, mm);
            if (mrw) {
                ms1 = tempf / sqrt(massi * massi);
                ms2 = tempf / sqrt(massj * massj);
                ms3 = -tempf / sqrt(massi * massj);
            }else{
                ms2 = ms1 = tempf;
                ms3 = -tempf;
            }
            basePoint[(eye + ll - 1) + (eye + mm - 1) * arrSz] += ms1;
            basePoint[(jay + ll - 1) + (jay + mm - 1) * arrSz] += ms2;
            basePoint[(eye + ll - 1) + (jay + mm - 1) * arrSz] += ms3;
            basePoint[(jay + ll - 1) + (eye + mm - 1) * arrSz] += ms3;
        }
    }
}

void Hessian::outputSparseMatrix(int i, int j, Real massi, Real massj, 
                                    Matrix3By3 rha, bool mrw, int arrSz, double *basePoint){

    for (int ll = 0; ll < 3; ll++)
        for (int mm = 0; mm < 3; mm++)
            if (mrw)
                basePoint[(i * 3 + ll) + (j * 3 + mm) * arrSz] +=
                        ((rha(ll, mm) / sqrt(massi * massj)));
            else 
                basePoint[(i * 3 + ll) + (j * 3 + mm) * arrSz] +=
                        rha(ll, mm);

}

void Hessian::evaluateCoulombBornRadiiPair(int i, int j, const Vector3DBlock *myPositions,
                                       const GenericTopology *myTopo,
                                       bool mrw, int mat_i, int mat_j, int mat_sz, double * mat_array) {
  int eye, jay;
  Real tempf, ms1, ms2, ms3;
  Matrix3By3 rha, rhc;
  ExclusionClass ec;

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
    //
    outputSparsePairMatrix(mat_i,mat_j,myTopo->atoms[i].scaledMass,myTopo->atoms[j].scaledMass,rha,mrw,mat_sz,mat_array);

  }

}

void Hessian::clear() {
  if (hessM != 0)
    for (unsigned int i = 0; i < sz * sz; i++) hessM[i] = 0.0;

}
}
