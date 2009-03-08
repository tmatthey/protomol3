#include <src/NormalModeProjection.h>

NormalModeProjection::NormalModeProjection() {
   mhQu = NULL;
}

NormalModeProjection::NormalModeProjection(int firstmode, int nummode, 
           Real gamma, int seed, Real temperature):
         StringNormalModeUtilities(firstmode, nummode, gamma, seed, temperature) {
    mhQu = NULL;
}

NormalModeProjection::~NormalModeProjection() {
   if (mhQu != NULL) delete [] mhQu;
}

void NormalModeProjection::initialize_utilities(int firstmode, 
          int nummode, Real gamma, int seed, Real temperature) {

    firstMode = firstmode;
    numMode = nummode;

    myGamma = gamma;

    mySeed = seed;

    myTemp = temperature;
}

void NormalModeProjection::CopyCommonBasis(double *ev, int b_d) {
   
    if (mhQu == NULL) mhQu = new double[_3N*_3N];
    _rfM = b_d;

    int _3N_2 = _3N*_rfM;

    for(int i=0;i<_3N_2;i++) mhQu[i] = ev[i];

    Q = &mhQu;

    cout<<"In NormalModeProjection :  _rfM "<<_rfM<<endl;

}

void NormalModeProjection::doModeSpaceProj(double *cc, Vector3DBlock *vi, Vector3DBlock *v0, int pts) {

   modeSpaceProj(cc,vi,v0);

   /*cout<<"NormalModeProjection ModeSpaceProj : Point "<<pts<<endl;
   for (int i=0;i<_rfM - (firstMode-1); i++) cout<<cc[i]<<" ";
   cout << endl; */

}
