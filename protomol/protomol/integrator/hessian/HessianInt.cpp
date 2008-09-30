#include <protomol/integrator/hessian/HessianInt.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/base/PMConstants.h>
#include <iostream>
#include <stdio.h>

#ifdef BUILD_FOR_FAH
#include <boost/iostreams/stream.hpp>
#include <fah/core/ChecksumDevice.h>

#else
#include <fstream>
#endif

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#else
#if defined (HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ HessianInt

const string HessianInt::keyword("hessianint");

HessianInt::HessianInt() :
  STSIntegrator() {
  eigVec = 0; eigVal = 0; eigIndx = 0;
}

HessianInt::HessianInt(Real timestep, string evec_s, string eval_s,
                       string hess_s, bool sorta, int fm, bool tef,
                       bool masswt, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), evecfile(evec_s),
  evalfile(eval_s), hessfile(hess_s), sortOnAbs(sorta), fixedModes(fm),
  textEig(tef), massWeight(masswt) {
  eigVec = 0; eigVal = 0; eigIndx = 0;
  //
  hsn.findForces(overloadedForces);         //find forces and parameters
}

HessianInt::~HessianInt() {
  int info;

  info = 0;
  //fix summed Hessian to average
  if (hsn.hessM != 0 && totStep > 1)
    for (unsigned int i = 0; i < sz * sz; i++)
      hsn.hessM[i] /= (double)totStep;

  //
  if (hsn.hessM != 0 && eigVec != 0 && eigVec != 0 && totStep &&
      (evecfile != "" || evalfile != "")) {

    report << hint << "[HessianInt::run] diagonalizing Hessian." << endr;
    //diagonalize
    info = diagHessian(eigVec, eigVal);
    if (info == 0) {
      int numneg;
      for (numneg = 0; numneg < (int)sz; numneg++)
        if (eigVal[numneg] > 0.0) break;

      //if(numneg) numneg--;
      report << hint << "[HessianInt::run] diagonalized! Number of negative "
        "eigenvalues = " << numneg << "." << endr;
      if (sortOnAbs) absSort();
      //output results
      //output eigenvec matrix/ eigenval vector
      outputDiagHess();
      report << hint << "[HessianInt::run] Nose Mass, Q = " << calcQ() << "."
             << endr;
    }else{
        report << hint << "[HessianInt::run] error info = " << info << endr;
    }
  }
  if (eigVec != 0) delete[] eigVec;
  if (eigVal != 0) delete[] eigVal;
  if (eigIndx != 0) delete[] eigIndx;
}

void HessianInt::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  //check options are sensible
  if (hessfile != "" && (evecfile != "" || evalfile != ""))
    THROW("[HessianInt::initialize] Cannot output Hessian after Lapack "
          "diagonalization!");

#if defined (HAVE_LAPACK)
#else
#if defined (HAVE_SIMTK_LAPACK)
#else
  if (evecfile != "" || evalfile != "")
    THROW("Hessian diagonalization requires Lapack libraries.");
#endif
#endif
  //creat hessian arrays and initialize sz
  sz = 3 * app->positions.size();
  report << hint << "[HessianInt::Find Hessian] sz=" << sz << endr;
  //
  hsn.initialData(sz);
  hsn.clear();
  //
  try{
    eigVec = new double[sz * sz];
    eigVal = new double[sz];
    //index for sort by absolute
    eigIndx = new int[sz];
  }catch(bad_alloc&){    
      report << error << "[HessianInt::initialize] Cannot allocate memory for "
             << "Eigenvectors." << endr;
  }
  for (unsigned int i = 0; i < sz; i++) eigIndx[i] = i;
  //
  totStep = 0;
}

void HessianInt::run(int numTimesteps) {
  if (numTimesteps < 1)
    return;

  if (totStep) {
    preStepModify();
    doHalfKickdoDrift();
    calculateForces();
    if (numTimesteps > 1)
      for (int i = 1; i < numTimesteps; i++) {
        doKickdoDrift();
        calculateForces();
        totStep++;
        //true for mass re-weight;
        hsn.evaluate(&app->positions, app->topology, massWeight);
      }

    else {
      totStep++;
      //true for mass re-weight;
      hsn.evaluate(&app->positions, app->topology, massWeight);
    }
    doHalfKick();
    postStepModify();
  } else {
    calculateForces();
    //Find current Hessian
    totStep++;
    //true for mass re-weight;
    hsn.evaluate(&app->positions, app->topology, massWeight);
    report << hint << "[HessianInt::Find Hessian] Hessian found!" << endr;
  }
  //
}

//Nose mass calculation based on Chris Sweet's Thesis
Real HessianInt::calcQ() {
    Real Q, sumF;
    
    if(sz > 6){
        sumF = 0;
        for(unsigned int i=6; i< sz;i++){
            if(eigVal[i]) sumF += sqrt(3.0 / fabs(eigVal[i]));
        }
        Q = 0.6/(8.0 * (sz - 6)) * sumF*sumF;
    }else{
        Q=0.0;
    }
    return Q;
}


void HessianInt::outputDiagHess() {
  unsigned int i;
#ifdef BUILD_FOR_FAH
    boost::iostreams::stream<FAH::ChecksumDevice> myFile;
#else
    ofstream myFile;
#endif

  if (hessfile != "") {
    myFile.open(hessfile.c_str(), ofstream::out);
    myFile.precision(10);
    //output hessian matrix to sparse form
    for (i = 0; i < sz * sz; i++)
      if (hsn.hessM[i] != 0.0)
        myFile << i / sz + 1 << " " << i % sz + 1
               << " " << hsn.hessM[i] << endl;

    //close file
    myFile.close();
  }
  if (evecfile != "") {
    //output eigenvec matrix
    myFile.open(evecfile.c_str(), ofstream::out);
    myFile.precision(10);
    //
    int vecnum = sz;
    int vecpos = sz / 3;
    //
    if (textEig) {
      myFile << vecnum * vecpos << endl;
      myFile << "! eigenvectors from Protomol/Lapack " << endl;
      for (int i = 0; i < vecnum; i++)
        for (int k = 0; k < vecpos; k++) {
          myFile << k + 1;
          for (int j = 0; j < 3; j++)
            myFile << " " << eigVec[eigIndx[i] * sz + k * 3 + j];

          myFile << endl;
        }

    } else {
      int numrec = sz / 3;
      int32 vp = vecpos;
      int32 fm = numrec * 3 - fixedModes;
      double ev = eigVal[eigIndx[(numrec - 1) * 3 + 2]];
      //
      //		myFile  << "! eigenvectors from Protomol/Lapack "<< endl;
      if (ISLITTLEENDIAN) swapBytes(vp);
      myFile.write((char *)&vp, sizeof(int32));
      if (ISLITTLEENDIAN) swapBytes(fm);
      myFile.write((char *)&fm, sizeof(int32));
      if (ISLITTLEENDIAN) swapBytes(ev);
      myFile.write((char *)&ev, sizeof(double));
      for (int i = 0; i < (int)(numrec * 3 - fixedModes); i++)
        for (int k = 0; k < vecpos; k++)
          //myFile << k + 1;
          for (int j = 0; j < 3; j++) {
            double evec = eigVec[eigIndx[i] * sz + k * 3 + j];
            if (ISLITTLEENDIAN) swapBytes(evec);
            myFile.write((char *)&evec, sizeof(double));
          }

    }
    //close file
    myFile.close();
  }
  if (evalfile != "") {
    //output eigenval vector
    myFile.open(evalfile.c_str(), ofstream::out);
    myFile.precision(10);
    //
    int numrec;
    numrec = sz / 3;
    myFile << numrec << endl;
    myFile << "! eigenvalues from Protomol/Lapack " << endl;
    for (int i = 0; i < numrec; i++) {
      myFile << i + 1;
      for (int j = 0; j < 3; j++)
        myFile << " " << eigVal[eigIndx[i * 3 + j]];

      myFile << endl;
    }

    //close file
    myFile.close();
  }
}

int HessianInt::diagHessian(double *eigVecO, double *eigValO) {   //
  double *wrkSp = 0;
  int *isuppz = 0, *iwork = 0;

  try{
    wrkSp = new double[26 * sz];
    isuppz = new int[2 * sz];
    iwork = new int[10 * sz];
  }catch(bad_alloc&){    
      report << error << "[HessianInt::diagHessian] Cannot allocate workspace "
             << "for diagonalization." << endr;
  }
  //Diagonalize
  int info = 0;                 /* output 0=success */
#if defined (HAVE_LAPACK) || defined (HAVE_SIMTK_LAPACK)
  /* LAPACK checks only first character N/V */
  char *jobz = "V"; char *range = "A"; char *uplo = "U";
  int n = sz;               /* order of coefficient matrix a  */
  int lda = sz;             /* leading dimension of a array*/
  double vl = 1.0;
  double vu = 1.0;
  int il = 1;
  int iu = 1;
  double abstol = 0;
  int m; int ldz = sz; int lwork = 26 * sz;   /* dimension of work array*/
  int liwork = 10 * sz;                       /* dimension of int work array*/
  //Recomended abstol for max precision
  char *cmach = "safe min";
#endif
  //call LAPACK
  //
#if defined ( HAVE_LAPACK )
  abstol = dlamch_(cmach);     //find machine safe minimum
  //
  dsyevr_(jobz, range, uplo, &n, hsn.hessM, &lda, &vl, &vu, &il, &iu, &abstol,
          &m, eigValO, eigVecO, &ldz, isuppz, wrkSp, &lwork, iwork, &liwork,
          &info);
#else
#if defined (HAVE_SIMTK_LAPACK)
  int len_cmach = 8;
  int len_jobz = 1; int len_range = 1; int len_uplo = 1;
  abstol = dlamch_(*cmach, len_cmach);     //find machine safe minimum
  //
  dsyevr_(*jobz, *range, *uplo, n, hsn.hessM, lda, &vl, &vu, &il, &iu,
          &abstol, m, eigValO, eigVecO, ldz, isuppz, wrkSp, lwork, iwork,
          &liwork, info, len_jobz, len_range, len_uplo);
#endif
#endif


  //delete arrays
  delete[] iwork;
  delete[] isuppz;
  delete[] wrkSp;
  //return status
  return info;
}

void HessianInt::doHalfKickdoDrift() {
  if (anyPreDriftOrNextModify()) {
    doHalfKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] += (*myForces)[i] * h * 0.5 /
                            app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void HessianInt::doKickdoDrift() {
  if (anyPreDriftOrNextModify() || anyPreStepModify() ||
      anyPostStepModify()) {
    if (anyPreStepModify() || anyPostStepModify()) {
      doHalfKick();
      postStepModify();
      preStepModify();
      doHalfKick();
    } else doKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void HessianInt::absSort() {    //sort for absolute magnitude
  unsigned int i;

  //find minimum abs value
  double minEv = fabs(eigVal[0]);
  for (i = 1; i < sz; i++)
    if (minEv < fabs(eigVal[i])) break;
    else minEv = fabs(eigVal[i]);

  i--;
  //sort around min
  if (i > 0) {
    int j = 0;
    eigIndx[j++] = i;
    int negp = i - 1;
    unsigned int posp = i + 1;
    while (negp >= 0 && posp < sz)
      if (fabs(eigVal[negp]) < fabs(eigVal[posp]))
        eigIndx[j++] = negp--;
      else eigIndx[j++] = posp++;

    while (negp >= 0) eigIndx[j++] = negp--;
  }
}

void HessianInt::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("eigvecFile",
               Value(evecfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Eigenvector filename")));
  parameters.push_back
    (Parameter("eigvalFile",
               Value(evalfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Eigenvalue filename")));
  parameters.push_back
    (Parameter("hessianFile",
               Value(hessfile, ConstraintValueType::NoConstraints()),
               string(""), Text("Hessian sparse filename")));
  parameters.push_back
    (Parameter("sortByAbs",
               Value(sortOnAbs, ConstraintValueType::NoConstraints()), false,
               Text("Sort eigenvalues/vectors by absolute magnitude")));
  parameters.push_back
    (Parameter("fixedModes",
               Value(fixedModes, ConstraintValueType::NoConstraints()),
               0, Text("Number of fixed modes")));
  parameters.push_back
    (Parameter("textEigFile",
               Value(textEig, ConstraintValueType::NoConstraints()), false,
               Text("Output eigenvectors as text file")));
  parameters.push_back
    (Parameter("massweight",
               Value(massWeight, ConstraintValueType::NoConstraints()),
               true, Text("Mass weight Hessian?")));
}

STSIntegrator *HessianInt::doMake(const vector<Value> &values,
                                  ForceGroup *fg) const {
  return new HessianInt(values[0], values[1], values[2], values[3], values[4],
                        values[5], values[6], values[7], fg);
}

