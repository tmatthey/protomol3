#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ NormalModeDiagonalize

const string NormalModeDiagonalize::keyword("NormalModeDiagonalize");

NormalModeDiagonalize::NormalModeDiagonalize() :
  MTSIntegrator(), NormalModeUtilities() {
  eigVec = 0; eigVal = 0; eigIndx = 0;
}

NormalModeDiagonalize::
NormalModeDiagonalize(int cycles, int avs, Real avss, int redi, int rayf,
                      int raya, int mins, Real minl,
                      ForceGroup *overloadedForces,
                      StandardIntegrator *nextIntegrator) :
  MTSIntegrator(cycles, overloadedForces, nextIntegrator),
  NormalModeUtilities(1, 1, 91.0, 1234, 300.0), noAvStep(avs),
  avStep(avss), rediagCount(redi), raylFrequ(rayf), raylAverage(raya),
  minSteps(mins), minLim(minl) {
  eigVec = 0; eigVal = 0; eigIndx = 0;
  //
  hsn.findForces(overloadedForces);         //find forces and parameters
  if (raylFrequ && raylFrequ < raylAverage) raylFrequ = raylAverage;
}

NormalModeDiagonalize::~NormalModeDiagonalize() {
  if (eigVec != 0 && eigAlloc) delete[] eigVec;
  if (eigVal != 0) delete[] eigVal;
  if (eigIndx != 0) delete[] eigIndx;
}

void NormalModeDiagonalize::initialize(ProtoMolApp *app) {
  MTSIntegrator::initialize(app);
  //test topmost integrator
  if (top() != this)
    report << error << "NormalModeDiagonalize not top integrator." << endr;
  //
  myNextNormalMode = dynamic_cast<NormalModeUtilities *>(myNextIntegrator);
  //Using complement of next integrator, so copy
  firstMode = myNextNormalMode->firstMode; numMode =
    myNextNormalMode->numMode;
  //NM initialization, but fix dof, as set by next integrator
  //last for complimentary forces, no gen noise
  NormalModeUtilities::initialize((int)app->positions.size(), app->topology,
                                  myForces, COMPLIMENT_FORCES);
  _rfM = myNextNormalMode->_rfM;
  //
  //do first force calculation, and remove non sub-space part
  //Need this or initial error, due to inner integrator energy?
  app->energies.clear();
  initializeForces();
  //
  //should check to see if mhQ is large enough and re-use, then fix
  //de-allocation
  if (mhQu != 0 && numEigvectsu == (unsigned int)_3N) {
    eigVec = mhQu;
    eigAlloc = false;
    validMaxEigv = true;
  } else {
    eigVec = new double[_3N * _3N];
    //mhQu = eigVec;
    eigAlloc = true;
    validMaxEigv = false;
  }
  eigVal = new double[_3N];
  eigIndx = new int[_3N];
  //initialize Hessian array
  hsn.initialData(_3N);
  if (eigVec == 0 || eigVal == 0 || eigIndx == 0 || hsn.hessM == 0)
    report << error << "Eigenvector array allocation error." << endr;
  //setup rediag/rayl counter incase valid
  nextRediag = rediagCount;
  nextRayl = raylFrequ;
  //clear do raylegh quo
  raylDo = false;
  //save positions where diagonalized for checkpoint save (assume I.C. if file)
  diagAt = *&app->positions;
}

//******************************************************************************
//****Normal run routine********************************************************
//******************************************************************************

void NormalModeDiagonalize::run(int numTimesteps) {
  if (numTimesteps < 1) return;

  if (rediagCount && raylDo == false &&
      (int)(app->topology->time / getTimestep()) >= nextRediag) {
    nextRediag += rediagCount;
    mhQu = 0;       //force rediag
  }
  //check valid eigenvectors
  if (mhQu == 0) {
    //Diagonalize if no input file
    report << hint <<
    "[NormalModeDiagonalize::run] Finding diagonalized Hessian." << endr;
    //save positions where diagonalized for checkpoint save
    diagAt = *&app->positions;
    Vector3DBlock tmpVel = *&app->velocities;
    //index for sort by absolute
    for (int i = 0; i < _3N; i++) eigIndx[i] = i;

    //mass re-weighted hessian
    hsn.clear();
    //true for mass re-weight;
    hsn.evaluate(&app->positions, app->topology, true);
    //minimized?
    if (minSteps > 0 && validMaxEigv) {
      //minimizer AND valid maximum eigenvalue
      mhQu = eigVec;        //use old eigs!
      minimizer(minLim, minSteps, true, false, true, &forceCalc, &lastLambda,
                &app->energies, &app->positions, app->topology);
      //SDminimize(minLim, minSteps,true);
      mhQu = 0;
    }
    //Avergaged?
    if (noAvStep > 1) {
      Real h = avStep * Constant::INV_TIMEFACTOR;
      const unsigned int count = app->positions.size();
      calculateForces();
      for (int nos = 1; nos < noAvStep; nos++) {
        report << debug(5) <<
        "[NormalModeDiagonalize::run] averaging step = " << nos <<
        " step size = " << avStep << endr;
        for (unsigned int i = 0; i < count; ++i) {
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
          app->positions[i] += app->velocities[i] * h;
        }

        //true for mass re-weight;
        hsn.evaluate(&app->positions, app->topology, true);
        calculateForces();
        for (unsigned int i = 0; i < count; ++i)
          app->velocities[i] += 
            (*myForces)[i] * 0.5 * h / app->topology->atoms[i].scaledMass;
      }

      for (int i = 0; i < _3N * _3N; i++) hsn.hessM[i] /= (double)noAvStep;

      //divide sum of Hessians
    }
    //average or minimized? then reset positions/velocities
    if (minSteps > 0 || noAvStep > 1) {
      *&app->positions = diagAt;        //back to original positions
      *&app->velocities = tmpVel;
    }
    //diagonalize
    int info = diagHessian(eigVec, eigVal, hsn.hessM);
    if (info == 0) {
      //find number of -ve eigs
      int ii;
      for (ii = 0; ii < _3N - 3; ii++) if (eigVal[ii + 3] > 0) break;

      //
      report 
        << hint << "[NormalModeDiagonalize::run] diagonalized. No. negative "
        "eigenvales = " << ii << endr;
      absSort(eigVec, eigVal, eigIndx);
    } else
      report << error << "Diagonalization failed." << endr;
    numEigvectsu = _3N;
    maxEigvalu = eigVal[_3N - 1];
    validMaxEigv = true;
    mhQu = eigVec;
    //sift current velocities/forces
    myNextNormalMode->subSpaceSift(&app->velocities, myForces);
    //If Rayleigh AND sifted, then calc QQ^T otherwise add eigvals to screen
    if (raylFrequ && eigAlloc) {
      report.precision(10);
      for (int ofs = 1; ofs < 5; ofs++) report << debug(4) <<
        "[NormalModeDiagonalize::run] EigVal " << ofs << " =" <<
        eigVal[_rfM - ofs] << endl;
    }
  }
  //do rayleigh quotient? Check allocated mhQ or too small for calc!
  if (raylFrequ && raylDo == false && eigAlloc == true &&
      (int)(app->topology->time / getTimestep()) >= nextRayl) {
    nextRayl += raylFrequ;
    raylDo = true;
    raylAvCount = 0;
    hsn.clear();        //clear Hessian
  }
  //main loop
  app->energies.clear();
  if (raylDo) {   //if calculating rayleigh do individulal steps
    int ii;
    for (ii = 0; ii < numTimesteps && raylAvCount < raylAverage;
         ii++, raylAvCount++) {
      myNextIntegrator->run(1);
      //true for mass re-weight;
      hsn.evaluate(&app->positions, app->topology, true);
    }

    if (ii < numTimesteps) myNextIntegrator->run(numTimesteps - ii);
    if (raylAvCount >= raylAverage) {
      //do calculation here, find fastest 5 quotients
      raylDo = false;
      double rQ, boundRq;
      for (int i = 0; i < 5; i++) {
        //calculate RQ from higest freq eig upwards
        calcRayleigh(&rQ, &boundRq, hsn.hessM, _rfM - (i + 1),
                     (double)raylAverage);
        
        report << debug(4) << "[NormalModeDiagonalize::run] rQ = " << rQ
               << " Bound=" << boundRq << " vector=" << _rfM - i << " _3N="
               << _3N << endl;
      }
    }
  } else
    myNextIntegrator->run(numTimesteps);
}

//******************************************************************************
//****Output int paramiters*****************************************************
//******************************************************************************

void NormalModeDiagonalize::getParameters(vector<Parameter> &parameters)
const {
  MTSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("averageSteps",
               Value(noAvStep, ConstraintValueType::NotNegative()), 1,
               Text("Hessian averaged over number of steps.")));
  parameters.push_back
    (Parameter("avStepSize", Value(avStep, ConstraintValueType::NotNegative()),
               1.0, Text("Step size for Hessian averaging.")));
  parameters.push_back
    (Parameter("reDiagFrequency",
               Value(rediagCount, ConstraintValueType::NotNegative()), 0,
               Text("Frequency of re-diagonalization (steps).")));
  parameters.push_back
    (Parameter("raylFrequency",
               Value(raylFrequ, ConstraintValueType::NotNegative()), 0,
               Text( "Frequency of Rayleigh Quationt calculation (steps).")));
  parameters.push_back
    (Parameter("raylAverage",
               Value(raylAverage, ConstraintValueType::NotNegative()), 1,
               Text("No. of steps to average Hessian for Rayleigh Quotient.")));
  parameters.push_back
    (Parameter("minSteps",
               Value(minSteps, ConstraintValueType::NotNegative()), 0,
               Text("Max. number of minimizer steps.")));
  parameters.push_back
    (Parameter("minLim", Value(minLim, ConstraintValueType::NotNegative()),
               1.0, Text("Minimization limit kcal mol^{-1}.")));
}

MTSIntegrator *NormalModeDiagonalize::
doMake(const vector<Value> &values, ForceGroup *fg,
       StandardIntegrator *nextIntegrator) const {
  return new NormalModeDiagonalize(values[0], values[1], values[2], values[3],
                                   values[4], values[5], values[6], values[7],
                                   fg, nextIntegrator);
}
