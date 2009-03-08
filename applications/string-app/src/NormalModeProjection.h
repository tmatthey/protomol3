/* -*- c++ -*- */
#ifndef NORMALMODEPROJECTION_H
#define NORMALMODEPROJECTION_H

#include <protomol/integrator/normal/NormalModeUtilities.h>

//
// Does projection from cartesian space to the subspace and back for all 
// intermediate points on the string. 
//

using namespace ProtoMol;


 class NormalModeProjection :  public NormalModeUtilities {

   public:
     NormalModeProjection();

     NormalModeProjection(int firstmode, int nummode, Real gamma, int seed, Real temperature);

     ~NormalModeProjection();

   public:

     void initialize_utilities(int firstmode, int nummode, Real gamma, int seed, Real temperature);

     void CopyCommonBasis(double *ev, int b_d);

     void CheckThisClass() {
       cout<<"NormalModeProjection : _rfM "<<_rfM<<", _3N "<<_3N<<", _N "<<_N<<", firstMode "<<firstMode<<endl;
     }

     void doModeSpaceProj(double *cc, Vector3DBlock *vi, Vector3DBlock *v0, int pts);


 };

#endif /* NORMALMODEPROJECTION_H */
