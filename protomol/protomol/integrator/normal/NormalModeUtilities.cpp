#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/ScalarStructure.h>

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#else
#if defined (HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

using namespace ProtoMol::Report;
using namespace ProtoMol;

//constructors
NormalModeUtilities::NormalModeUtilities() :
  firstMode(1), numMode(-1), myGamma(-1), mySeed(-1), myTemp(-1) {
  tmpFX = 0; tmpC = 0; invSqrtMass = 0; sqrtMass = 0;
}

NormalModeUtilities::NormalModeUtilities(int firstmode, int nummode,
                                         Real gamma, int seed,
                                         Real temperature) :
  firstMode(firstmode), numMode(nummode),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)), mySeed(seed),
  myTemp(temperature) {
  tmpFX = 0; tmpC = 0; invSqrtMass = 0; sqrtMass = 0;
}

NormalModeUtilities::~NormalModeUtilities() {
  //
  if (tmpFX != 0) delete[] tmpFX;
  if (tmpC != 0) delete[] tmpC;
  if (invSqrtMass != 0) delete[] invSqrtMass;
  if (sqrtMass != 0) delete[] sqrtMass;
}

void NormalModeUtilities::initialize(int sz, GenericTopology *myTopo,
                                     Vector3DBlock *myForces,
                                     int nm_flags) {
#if defined (HAVE_LAPACK)
#else
#if defined (HAVE_SIMTK_LAPACK)
#else
  report << error << "Normal Mode integrators require Lapack libraries." 
         << endr;
#endif
#endif
  //local force pointer
  myForcesP = myForces;
  //type
  if (nm_flags & COMPLIMENT_FORCES) complimentForces = true;
  else complimentForces = false;
  if (nm_flags & GEN_COMP_NOISE) genCompNoise = true;
  else genCompNoise = false;
  //
  _N = sz;
  _3N = 3 * _N;   //used everywhere, must be initialized
  //setup for auto re-diag by allowing full set of eigenvectors
  if (!numEigvectsu) numEigvectsu = _3N;
  //first mode?
  if (firstMode < 1 || firstMode > _3N)
    report << error << "firstmode = " << firstMode
           << ", must be > 0 and less than 3N = " << _3N << "." << endr;
  //check numMode
  if (numMode < 1 || numMode > _3N - (firstMode - 1))
    report 
      << error << "numbermodes = " << numMode 
      << ", must be > 0 and less than or equal to 3N-(firstmode-1) = "
      << _3N - (firstMode - 1) << "." << endr;
  //
  if (numMode > (int)numEigvectsu - (firstMode - 1))
    report 
      << error << "Insufficient eigenvectors: numbermodes = "
      << numMode << ", must be less than than or equal to m-(firstmode-1) = "
      << numEigvectsu - (firstMode - 1) << "." << endr;

  //number of low frequency modes
  _rfM = numMode + (firstMode - 1);
  //Degrees of freedom
  myTopo->degreesOfFreedom = _rfM - (firstMode - 1);
  if (firstMode == 1) myTopo->degreesOfFreedom -= 3;
  //temporary mode space variable for intermediate calculations
  tmpC = new double[_3N];
  //define temporary position/force vector
  tmpFX = new double[_3N];
  //setup sqrt masses and their inverse
  invSqrtMass = new double[_N];
  sqrtMass = new double[_N];
  for (int i = 0; i < _N; i++) {
    sqrtMass[i] = sqrt(myTopo->atoms[i].scaledMass);
    if (sqrtMass[i]) invSqrtMass[i] = 1.0 / sqrtMass[i];
  }

  //
  posTemp.resize(_N);
  //***********************************************************
  //gauss randoms
  gaussRandCoord1.resize(_N);
  gaussRandCoord2.resize(_N);
  if (genCompNoise) tempV3DBlk.resize(_N);  //use comliment noise gen if flagged
}

//test for remaining modes, equivalent to old fixModes=0
bool NormalModeUtilities::testRemainingModes() {
  if ((numMode + (firstMode - 1)) < _3N) return true;
  else return false;
}

//Setup string of integrators so that all point to the top integrator
// 'Eigenpointers' :-)
void NormalModeUtilities::setIntegratorSetPointers(Integrator *integrator,
                                                   EigenvectorInfo *eipt,
                                                   bool eiValid) {
  NormalModeUtilities *nmint;
  if (eiValid) {
    mhQu = eipt->myEigenvectors;
    maxEigvalu = eipt->myMaxEigenvalue;
    numEigvectsu = eipt->myNumEigenvectors;
  } else {
    mhQu = 0;
    maxEigvalu = 0;
    numEigvectsu = 0;
  }
  eigValP = &maxEigvalu;
  Q = &mhQu;
  //find next integrators
  for (Integrator *i = integrator->next(); i != 0; i = i->next()) {
    nmint = dynamic_cast<NormalModeUtilities *>(i);
    if (nmint == 0) report << error <<
      "Normal Mode integrator chain contains unknown integrator type." <<
      endr;
    nmint->eigValP = &maxEigvalu;
    nmint->Q = &mhQu;
    nmint->numEigvectsu = numEigvectsu;
    nmint->mhQu = 0;
    nmint->maxEigvalu = 0;
  }
}

//****Routines to convert between Vector3DBlocks and linear arrays**************

//Convert Vector3DBlock formal to linear array for BLAS
double *NormalModeUtilities::vector3DBlockTOvect(Vector3DBlock *blkDat,
                                                 double *vecDat) {
  for (int i = 0; i < 3 * _N; i++) vecDat[i] = (*blkDat)[i / 3][i % 3];

  return vecDat;
}

//Convert linear array from BLAS to Vector3DBlock
Vector3DBlock *NormalModeUtilities::vectTOvector3DBlock(double *vecDat,
                                                        Vector3DBlock *blkDat)
{
  for (int i = 0; i < 3 * _N; i++) (*blkDat)[i / 3][i % 3] = vecDat[i];

  return blkDat;
}

//******************************************************************************
//****Projectors for complement sub space and sub space*************************
//******************************************************************************

//Find forces acting outside subspace
Vector3DBlock *NormalModeUtilities::nonSubspaceForce(Vector3DBlock *force,
                                                     Vector3DBlock *iPforce) {
  //
  vector3DBlockTOvect(force, tmpFX);    //put positions-x_0 into linear array
  //calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
  //f'=M^{-1/2}*f
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= invSqrtMass[i / 3];

  //c=hQ^T*M^{-1/2}*f
  double alpha = 1.0; double beta = 0.0;        //multiplyers, see Blas docs.
#if defined (HAVE_LAPACK)
  char *transA = "T";     // Transpose Q, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM; int incxy = 1;     //sizes
  dgemv_(transA, &m, &n, &alpha, (*Q), &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transA = "T";     // Transpose Q, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM; int incxy = 1;     //sizes
  int len_transa = 1;                           //length of transA
  dgemv_(*transA, m, n, alpha, (*Q), m, tmpFX, incxy, beta, tmpC, incxy,
         len_transa);
#endif
#endif
  //
  //calculate f''=M^{-1/2}*f'-hQc using BLAS
  alpha = -1.0;   beta = 1.0;
  //
#if defined (HAVE_LAPACK)
  char *transB = "N";
  dgemv_(transB, &m, &n, &alpha, (*Q), &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transB = "N";
  dgemv_(*transB, m, n, alpha, (*Q), m, tmpC, incxy, beta, tmpFX, incxy,
         len_transa);
#endif
#endif
  //
  //f'''=M^{1/2}*f''
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= sqrtMass[i / 3];

  //put back into vector3DBlocks
  vectTOvector3DBlock(tmpFX, iPforce);
  //delete temporary array
  return iPforce;
}

//Find positions outside subspace
Vector3DBlock *NormalModeUtilities::nonSubspacePosition(Vector3DBlock *force,
                                                        Vector3DBlock *iPforce)
{
  //
  vector3DBlockTOvect(force, tmpFX);    //put positions-x_0 into linear array
  //calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
  //f'=M^{-1/2}*f
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= sqrtMass[i / 3];

  //c=hQ^T*M^{-1/2}*f
  double alpha = 1.0; double beta = 0.0;        //multiplyers, see Blas docs.
  //
#if defined (HAVE_LAPACK)
  char *transA = "T";     // Transpose Q, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM; int incxy = 1;     //sizes
  dgemv_(transA, &m, &n, &alpha, (*Q), &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transA = "T";     // Transpose Q, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM; int incxy = 1;     //sizes
  int len_transa = 1;                           //length of transA
  dgemv_(*transA, m, n, alpha, (*Q), m, tmpFX, incxy, beta, tmpC, incxy,
         len_transa);
#endif
#endif
  //
  //calculate f''=M^{-1/2}*f'-hQc using BLAS
  alpha = -1.0;   beta = 1.0;
  //
#if defined (HAVE_LAPACK)
  char *transB = "N";
  dgemv_(transB, &m, &n, &alpha, (*Q), &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transB = "N";
  dgemv_(*transB, m, n, alpha, (*Q), m, tmpC, incxy, beta, tmpFX, incxy,
         len_transa);
#endif
#endif
  //
  //f'''=M^{1/2}*f''
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= invSqrtMass[i / 3];

  //put back into vector3DBlocks
  vectTOvector3DBlock(tmpFX, iPforce);
  //delete temporary array
  return iPforce;
}

//Find forces acting inside subspace
Vector3DBlock *NormalModeUtilities::subspaceForce(Vector3DBlock *force,
                                                  Vector3DBlock *iPforce) {
  //
  vector3DBlockTOvect(force, tmpFX);    //put positions-x_0 into linear array
  //calculate M^{1/2}QQ^TM^{-1/2}f using BLAS
  //f'=M^{-1/2}*f
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= invSqrtMass[i / 3];

  //c=Q^T*M^{-1/2}*f
  double alpha = 1.0; double beta = 0.0;
  //
#if defined (HAVE_LAPACK)
  char *transA = "T";      // Transpose, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM - (firstMode - 1); int incxy = 1;   //sizes
  dgemv_(transA, &m, &n, &alpha, &((*Q)[_3N * (firstMode - 1)]), &m, tmpFX,
         &incxy, &beta, tmpC,
         &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transA = "T";      // Transpose, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM - (firstMode - 1); int incxy = 1;   //sizes
  int len_transa = 1;
  dgemv_(*transA, m, n, alpha, &((*Q)[_3N * (firstMode - 1)]), m, tmpFX,
         incxy, beta, tmpC, incxy,
         len_transa);
#endif
#endif
  //
  //f''=Qc
  alpha = 1.0;    beta = 0.0;
  //
#if defined (HAVE_LAPACK)
  char *transB = "N";  /* LAPACK checks only first character N/V */
  dgemv_(transB, &m, &n, &alpha, &((*Q)[_3N * (firstMode - 1)]), &m, tmpC,
         &incxy, &beta, tmpFX,
         &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transB = "N";  /* LAPACK checks only first character N/V */
  dgemv_(*transB, m, n, alpha, &((*Q)[_3N * (firstMode - 1)]), m, tmpC, incxy,
         beta, tmpFX, incxy,
         len_transa);
#endif
#endif
  //f'''=M^{1/2}*f''
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= sqrtMass[i / 3];

  //put back into vector3DBlocks
  vectTOvector3DBlock(tmpFX, iPforce);
  //delete temporary array
  return iPforce;
}

//Find velocities acting inside subspace
Vector3DBlock *NormalModeUtilities::subspaceVelocity(Vector3DBlock *force,
                                                     Vector3DBlock *iPforce) {
  //
  vector3DBlockTOvect(force, tmpFX);    //put positions-x_0 into linear array
  //calculate M^{-1/2}QQ^TM^{1/2}f using BLAS
  //v'=M^{1/2}*v
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= sqrtMass[i / 3];

  //c=Q^T*M^{-1/2}*v
  double alpha = 1.0; double beta = 0.0;
  //
#if defined (HAVE_LAPACK)
  char *transA = "T";        //Transpose, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM - (firstMode - 1); int incxy = 1;   //sizes
  dgemv_(transA, &m, &n, &alpha, &((*Q)[_3N * (firstMode - 1)]), &m, tmpFX,
         &incxy, &beta, tmpC,
         &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transA = "T";        //Transpose, LAPACK checks only first character N/V
  int m = _3N; int n = _rfM - (firstMode - 1); int incxy = 1;   //sizes
  int len_transa = 1;
  dgemv_(*transA, m, n, alpha, &((*Q)[_3N * (firstMode - 1)]), m, tmpFX,
         incxy, beta, tmpC, incxy,
         len_transa);
#endif
#endif
  //
  //v''=Qc
  alpha = 1.0;    beta = 0.0;
  //
#if defined (HAVE_LAPACK)
  char *transB = "N";  /* LAPACK checks only first character N/V */
  dgemv_(transB, &m, &n, &alpha, &((*Q)[_3N * (firstMode - 1)]), &m, tmpC,
         &incxy, &beta, tmpFX,
         &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  char *transB = "N";  /* LAPACK checks only first character N/V */
  dgemv_(*transB, m, n, alpha, &((*Q)[_3N * (firstMode - 1)]), m, tmpC, incxy,
         beta, tmpFX, incxy,
         len_transa);
#endif
#endif
  //v'''=M^{-1/2}*v''
  for (int i = 0; i < _3N; i++)
    tmpFX[i] *= invSqrtMass[i / 3];

  //put back into vector3DBlocks
  vectTOvector3DBlock(tmpFX, iPforce);
  //delete temporary array
  return iPforce;
}

void NormalModeUtilities::subSpaceSift(Vector3DBlock *velocities,
                                       Vector3DBlock *forces) {
  //sift current data into subspace
  subspaceVelocity(velocities, velocities);
  subspaceForce(forces, forces);
}

//force projection for post force modifier
void NormalModeUtilities::forceProjection() {
  if ((*Q) != 0)
    if (complimentForces) nonSubspaceForce(myForcesP, myForcesP);
    else subspaceForce(myForcesP, myForcesP);
}

//******************************************************************************
//****Langevin Thermostat*******************************************************
//******************************************************************************

// Generate projected vector of gausians
void NormalModeUtilities::genProjGauss(Vector3DBlock *gaussRandCoord,
                                       GenericTopology *myTopo) {
  //generate set of random force variables and project into sub space
  if ((int)gaussRandCoord->size() != _N) return;
  for (int i = 0; i < _3N; i++)
    (*gaussRandCoord)[i / 3][i % 3] = randomGaussianNumber(mySeed);

  //
  for (int i = 0; i < _N; i++)
    (*gaussRandCoord)[i] *= sqrtMass[i];

  if (complimentForces) nonSubspaceForce(gaussRandCoord, gaussRandCoord);
  else subspaceForce(gaussRandCoord, gaussRandCoord);
  for (int i = 0; i < _N; i++)
    (*gaussRandCoord)[i] /= myTopo->atoms[i].scaledMass;
}

// Generate projected vector of gausians AND its compliment
void NormalModeUtilities::genProjGaussC(Vector3DBlock *gaussRandCoord,
                                        Vector3DBlock *gaussRandCoordm,
                                        GenericTopology *myTopo) {
  //generate set of random force variables and project into sub space
  if ((int)gaussRandCoord->size() != _N || (int)gaussRandCoordm->size() !=
      _N) return;
  for (int i = 0; i < _3N; i++)
    (*gaussRandCoord)[i / 3][i % 3] = randomGaussianNumber(mySeed);

  //
  for (int i = 0; i < _N; i++)
    (*gaussRandCoord)[i] *= sqrtMass[i];

  //get randoms for compliment
  *gaussRandCoordm = *gaussRandCoord;
  //
  if (complimentForces) nonSubspaceForce(gaussRandCoord, gaussRandCoord);
  else subspaceForce(gaussRandCoord, gaussRandCoord);
  //
  for (int i = 0; i < _N; i++) {
    (*gaussRandCoordm)[i] -= (*gaussRandCoord)[i];
    (*gaussRandCoord)[i] /= myTopo->atoms[i].scaledMass;
    (*gaussRandCoordm)[i] /= myTopo->atoms[i].scaledMass;
  }
}

// fluctuation using Dr. Skeel's LI scheme which involves a semi-update
// of velocities and a complete update of positions
void NormalModeUtilities::nmlDrift(Vector3DBlock *myPositions,
                                   Vector3DBlock *myVelocities, Real dt,
                                   GenericTopology *myTopo) {
  //const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
  const Real tau1 = (1.0 - exp(-myGamma * dt)) / myGamma;
  const Real tau2 = (1.0 - exp(-2 * myGamma * dt)) / (2 * myGamma);
  const Real forceConstant = 2 * Constant::BOLTZMANN * myTemp * myGamma;
  const Real sqrtTau2 = sqrt(tau2);

  //  It is possible that roundoff and/or truncation can make this value < 0.
  //  It should be set to 0 if this happens.
  Real sqrtVal1 = (dt - tau1 * tau1 / tau2);

  if (sqrtVal1 < 0.)
    sqrtVal1 = 0;
  else
    sqrtVal1 = sqrt(sqrtVal1);

  Real sqrtFCoverM = sqrt(forceConstant);   /// mass );
  Real langDriftVal = sqrtFCoverM / myGamma;
  Real langDriftZ1 = langDriftVal * (tau1 - tau2) / sqrtTau2;
  Real langDriftZ2 = langDriftVal * sqrtVal1;

  //generate 1st and 2nd set of random force variables and project into sub 
  //space and comliment if required
  if (genCompNoise) genProjGaussC(&gaussRandCoord1, &tempV3DBlk, myTopo);
  else genProjGauss(&gaussRandCoord1, myTopo);
  genProjGauss(&gaussRandCoord2, myTopo);
  // Update positions and correct (semi-update) velocity
  for (unsigned int i = 0; i < myPositions->size(); i++) {
    posTemp[i] =
      (gaussRandCoord1[i] * langDriftZ1 + gaussRandCoord2[i] * langDriftZ2 +
       (*myVelocities)[i]) * tau1;
    // semi-update velocities
    (*myVelocities)[i] = (*myVelocities)[i] * exp(-myGamma * dt) +
      gaussRandCoord1[i] * sqrtFCoverM * sqrtTau2;
  }

  //Check position change is in the sub space and add to positions
  subspaceVelocity(&posTemp, &posTemp);
  for (unsigned int i = 0; i < myPositions->size(); i++)
    (*myPositions)[i] += posTemp[i];

  //fix COM/momentum (not conserved)
  buildMolecularCenterOfMass(myPositions, myTopo);
  buildMolecularMomentum(myVelocities, myTopo);
}

//******************************************************************************
//****Diagonalization/ minimizer etc.*******************************************
//******************************************************************************

//Diagonalize hessian***********************************************************
int NormalModeUtilities::diagHessian(double *eigVecO, double *eigValO,
                                     double *hsnhessM) {
  double *wrkSp;
  int *isuppz, *iwork;

  wrkSp = new double[26 * _3N];
  isuppz = new int[2 * _3N];
  iwork = new int[10 * _3N];
  //Diagonalize
  /* LAPACK checks only first character N/V */
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
  char *jobz = "V"; char *range = "A"; char *uplo = "U"; 
  int n = _3N;               /* order of coefficient matrix a  */
  int lda = _3N;             /* leading dimension of a array*/
  double vl = 1.0;
  double vu = 1.0;
  int il = 1;
  int iu = 1;
  double abstol = 0;
  int m; int ldz = _3N; int lwork = 26 * _3N;   /* dimension of work array*/
  int liwork = 10 * _3N;                        /* dimension of int work array*/
  //Recomended abstol for max precision
  char *cmach = "safe min";
#endif
  int info = 0;                 /* output 0=success */
  //call LAPACK
  //
#if defined ( HAVE_LAPACK )
  abstol = dlamch_(cmach);     //find machine safe minimum
  //
  dsyevr_(jobz, range, uplo, &n, hsnhessM, &lda, &vl, &vu, &il, &iu, &abstol,
          &m, eigValO, eigVecO, &ldz, isuppz,
          wrkSp, &lwork, iwork, &liwork,
          &info);
#else
#if defined (HAVE_SIMTK_LAPACK)
  int len_cmach = 8;
  int len_jobz = 1; int len_range = 1; int len_uplo = 1;
  abstol = dlamch_(*cmach, len_cmach);     //find machine safe minimum
  //
  dsyevr_(*jobz, *range, *uplo, n, hsnhessM, lda, &vl, &vu, &il, &iu, &abstol,
          m, eigValO, eigVecO, ldz, isuppz,
          wrkSp, lwork, iwork, &liwork, info, len_jobz, len_range,
          len_uplo);
#endif
#endif
  //delete arrays
  delete[] iwork;
  delete[] isuppz;
  delete[] wrkSp;
  //return status
  return info;
}

//sort vectors for absolute value*********************************************//
void NormalModeUtilities::absSort(double *eigVec, double *eigVal,
                                  int *eigIndx) {
  int i;

  //find minimum abs value
  double minEv = fabs(eigVal[0]);
  for (i = 1; i < _3N; i++)
    if (minEv < fabs(eigVal[i])) break;
    else minEv = fabs(eigVal[i]);

  i--;
  //sort around min
  if (i > 0) {
    int j = 0;
    eigIndx[j++] = i;
    int negp = i - 1;
    int posp = i + 1;
    while (negp >= 0 && posp < _3N)
      if (fabs(eigVal[negp]) < fabs(eigVal[posp]))
        eigIndx[j++] = negp--;
      else eigIndx[j++] = posp++;

    while (negp >= 0) eigIndx[j++] = negp--;
  }
  //Sort actual eigenvector array
  double *tmpVect = new double[_3N];
  double tmpElt;
  int ii, k;
  for (i = 0; i < _3N; i++)
    if (eigIndx[i] != (int)i && eigIndx[i] != -1) {         //need to swap?
      for (int j = 0; j < _3N; j++) {
        tmpVect[j] = eigVec[i * _3N + j];
        eigVec[i * _3N + j] = eigVec[eigIndx[i] * _3N + j];
      }

      eigIndx[i] = -1;                                      //flag swapped
      ii = i;
      do {
        for (k = 0; k < _3N && eigIndx[k] != (int)ii; k++) ;

        //find where tmpVect goes
        if (k == _3N || k == ii) break;     //end chain where indeces are equal
        for (int j = 0; j < _3N; j++) {             //put it there
          tmpElt = tmpVect[j];
          tmpVect[j] = eigVec[k * _3N + j];
          eigVec[k * _3N + j] = tmpElt;
        }

        eigIndx[k] = -1;                                    //flag swapped
        ii = k;
      } while (k < _3N);
    }

  delete[] tmpVect;
}

//calculate the Rayleigh quotient
double NormalModeUtilities::calcRayleigh(double *rQ, double *boundRq,
                                         double *hsnhessM, int numv,
                                         double raylAverage) {
  //Norm Mode diagnostics-Rayleigh Quotient
  //Get current Hessian
  for (int i = 0; i < _3N * _3N; i++) hsnhessM[i] /= raylAverage;

  //divide sum of Hessians
  //set BLAS variables
  //calculate Rayleigh Quo/bound
  //double rQ, boundRq;
  //define temporary position/force vector
  //
#if defined (HAVE_LAPACK)
  double alpha = 1.0; double beta = 0.0;
  char *transA = "N";
  int m = _3N; int n = _3N; int incxy = 1;    //sizes
  double *tmpFX = new double[_3N];
  dgemv_(transA, &m, &n, &alpha, hsnhessM, &m, &((*Q)[_3N * (numv)]), &incxy,
         &beta, tmpFX,
         &incxy);
  *rQ = ddot_(&n, &((*Q)[_3N * (numv)]), &incxy, tmpFX, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  double alpha = 1.0; double beta = 0.0;
  char *transA = "N";
  int m = _3N; int n = _3N; int incxy = 1;    //sizes
  double *tmpFX = new double[_3N];
  int len_transa = 1;
  dgemv_(*transA, m, n, alpha, hsnhessM, m, &((*Q)[_3N * (numv)]), incxy,
         beta, tmpFX, incxy,
         len_transa);
  *rQ = ddot_(n, &((*Q)[_3N * (numv)]), incxy, tmpFX, incxy);
#endif
#endif
  //bound
  for (int i = 0; i < _3N; i++) hsnhessM[i + _3N * i] -= *rQ;

#if defined (HAVE_LAPACK)
  dgemv_(transA, &m, &n, &alpha, hsnhessM, &m, &((*Q)[_3N * (numv)]), &incxy,
         &beta, tmpFX,
         &incxy);
  *boundRq = dnrm2_(&n, tmpFX, &incxy);
#else
#if defined (HAVE_SIMTK_LAPACK)
  dgemv_(*transA, m, n, alpha, hsnhessM, m, &((*Q)[_3N * (numv)]), incxy,
         beta, tmpFX, incxy,
         len_transa);
  *boundRq = dnrm2_(n, tmpFX, incxy);
#endif
#endif
  for (int i = 0; i < _3N; i++) hsnhessM[i + _3N * i] += *rQ;

  //restore
  //
  delete[] tmpFX;
  //
  return *rQ;
}

//******************************************************************************
//****Minimizer*****************************************************************
//******************************************************************************

//Simple-Steepest decent minimizer for all modes outside subspace.
//PK update respects 'c' subspace positions. Requires virtual force calculation
//function utilityCalculateForces().
int NormalModeUtilities::minimizer(Real peLim, int numloop, bool simpM,
                                   bool reDiag, bool nonSubspace,
                                   int *forceCalc, Real *lastLambda,
                                   ScalarStructure *myEnergies,
                                   Vector3DBlock *myPositions,
                                   GenericTopology *myTopo) {
  int in, itr, numLambda;
  Real oldPot, lambda, lambda1, lambdaSlp, lambdaSlp1, lastDiff;
  int rsCG;

  //initialize
  rsCG = 0; *lastLambda = 0.0; numLambda = 0; itr = 0; *forceCalc = 0;
  lastDiff = 5; //start at 5 kcal mol^{-1}
  //find new forces
  utilityCalculateForces();
  (*forceCalc)++;
  //Set first value of /lambda to be 1/eigval
  //exact solution for highest frequency mode if force mass weighted
  lambda = 1.0 / *eigValP;
  for (in = 0; in < numloop; in++) {
    itr++;
    //
    report.precision(10);
    report << debug(6) << "[NormalModeUtilities::minimizer] PE= " <<
    myEnergies->potentialEnergy() << endr;
    //****find search direction vector posTemp
    //find forces in compliment space.
    if (nonSubspace) nonSubspaceForce(myForcesP, myForcesP);
    //sift so position move is in compliment space, mass weighted. Replaces 
    //nonSubspacePosition(myForces, myForces)
    for (int i = 0; i < _N; i++) (*myForcesP)[i] /=
        myTopo->atoms[i].scaledMass;

    //set posTemp
    posTemp = *myForcesP;
    //find slope of original PE with /lambda=0 here.
    lambdaSlp1 = 0.0;
    for (int k = 0; k < _N; k++) lambdaSlp1 -= posTemp[k].dot((*myForcesP)[k]);

    //save PE at /lambda=0
    oldPot = myEnergies->potentialEnergy();
    report << debug(7) << "[NormalModeUtilities::minimizer] lambd= " <<
    lambda << endl;
    //find force at new position pos+/lambda*posTemp
    (*myPositions).intoWeightedAdd(lambda, posTemp);
    utilityCalculateForces();
    (*forceCalc)++;
    //Full minimizer? then solve for quadratic minimum
    if (!simpM) {
      //find slope of PE with /lambda here
      lambdaSlp = 0.0;
      for (int k = 0; k < _N; k++) lambdaSlp -= posTemp[k].dot(
          (*myForcesP)[k]);

      //solve for minimum for quadratic fit using two PE vales and the slope 
      //with /lambda
      Real a, b, oldLambda;
      oldLambda = lambda;
      a =
        -((myEnergies->potentialEnergy() -
           oldPot) / lambda - lambdaSlp) / lambda;
      b = lambdaSlp - 2.0 * a * lambda;
      lambda = -b / (2 * a);
      if (lambda <= 0.0) lambda = oldLambda;
      else {
        //Put solution into positions (but remove temporary solution for 
        //quadratic fit via oldLambda)
        (*myPositions).intoWeightedAdd(lambda - oldLambda, posTemp);
        utilityCalculateForces();
        (*forceCalc)++;
      }
    }
    //end full minimizer
    //update total gamma
    *lastLambda += lambda;
    numLambda++;
    //test for end, too large lambda test first
    if ((oldPot - myEnergies->potentialEnergy()) < 0)
      if (rsCG > 4)
        report
          << error
          << "[NormalModeUtilities::minimizer] Minimization failed, Aborting." 
          << endr;
      else
      if (!reDiag) {          //allow minimization if mode at angle to sub-space
        //calc optimum lambda from first slope
        Real a1;
        a1 = (myEnergies->potentialEnergy() - oldPot - lambdaSlp1 * lambda) /
          (lambda * lambda);
        lambda1 = -lambdaSlp1 / (2 * a1);
        //Test that the quadratic solution gives a predicted PE value
        //where the difference from the old PE value is bounded by the 
        //difference from the last succesful step, else solve quadratic for
        //the last difference.
        Real calcPE = 
          a1 * lambda1 * lambda1 + lambdaSlp1 * lambda1 + oldPot;
        if (oldPot - calcPE > lastDiff && lambdaSlp1 != 0.0)
          lambda1 = (-lastDiff * 2.0) / lambdaSlp1;
        (*myPositions).intoWeightedAdd(-lambda, posTemp);    //reset positions
        *lastLambda -= lambda;
        numLambda--;
        utilityCalculateForces();
        (*forceCalc)++;
        if (lambda1 > 0.0 && lambda1 < lambda) lambda = lambda1;
        else lambda /= 2.0;
        rsCG++;
        report 
          << debug(1)
          << "[NormalModeUtilities::minimizer] Reset CG, PE fail. Cycle= "
          << rsCG << " lambda= " << lambda << endl;
      } else {
        (*myPositions).intoWeightedAdd(-lambda, posTemp);      //reset positions
        (*Q) = 0;               //force rediagonalization
        return -1;                          //flag aborted
      }
    else {
      rsCG = 0;
      lambda = 1.0 / *eigValP;          //revert to original value
    }
    if ((oldPot - myEnergies->potentialEnergy()) < peLim && !rsCG) break;
    if (!rsCG) lastDiff = oldPot - myEnergies->potentialEnergy();
    //
  }

  if (numLambda) *lastLambda /= (Real)numLambda;
  else lastLambda = 0;
  return itr;
}

