/* -*- c++ -*- */
#ifndef HESSIAN_H
#define HESSIAN_H

#include <protomol/force/Force.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  /**
   *
   * Calculates the Hessian or mass re-weighted hessian
   * for the current force field.
   *
   */
  class Hessian {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Hessian();
    Hessian(unsigned int szin);   //dimension of matrix (sz*sz)
    Hessian(const Hessian &hess); //copy const.
    ~Hessian();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Hessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialData(unsigned int szin);
    void findForces(ForceGroup *overloadedForces);
    void evaluate(const Vector3DBlock *myPositions,   //positions
                  const GenericTopology *myTopo,      //topology
                  const bool mrw);                    //mass re-weighted

    void evaluateCoulombSCPISM(const Vector3DBlock *myPositions,
                               const GenericTopology *myTopo,
                               bool mrw);
    void evaluateCoulombBornRadii(const Vector3DBlock *myPositions,
                                  const GenericTopology *myTopo,
                                  bool mrw);

  public:
    void clear(); // clear the hessian matrix

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class Hessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    bool myBond;    //force field components
    bool myAngle;
    bool myCoulomb;
    bool myCoulombDielec;
    bool myCoulombSCPISM;
    bool myCoulombBornRadii;
    bool myLennardJones;
    bool myDihedral;
    bool myImproper;
    Real cCutoff, cSwitchon, cSwitch, cOrder, cSwitchoff;   //switch data
    Real lCutoff, lSwitchon, lSwitch, lOrder, lSwitchoff;
    Real D, S, epsi;
    unsigned int sz; //size
    int swt;
  public:
    double *hessM;  //matrix
  };
}
#endif
