/* -*- c++ -*- */
#ifndef COULOMBBORNRADIIFORCE_H
#define COULOMBBORNRADIIFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/type/Real.h>
#include <protomol/type/Vector3D.h>
#include <protomol/config/Parameter.h>
#include <string>

namespace ProtoMol {
  struct Internal {
    Internal(int i, Vector3D v, Real r) : atom2(i), diff(v), dist(r) {}
    
    int atom2;
    Vector3D diff;
    Real dist;
  };
  
  struct InteractionInfo {
    vector<Internal> internal;
    vector<Real> dRs;
  };
  
  static vector<InteractionInfo> info;

  //____ CoulombBornRadiiForce
  class CoulombBornRadiiForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom, so the Coulomb
    // constant here is one.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CoulombBornRadiiForce() {};
    CoulombBornRadiiForce(int s) {sw = s;};

  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombBornRadiiForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) {
      sw = topo->doSCPISM;
      Real q_i = topo->atoms[atom1].scaledCharge;
      Real q_j = topo->atoms[atom2].scaledCharge;

      if ((q_i == 0) || (q_j == 0)) return;

      energy = 0.0;
      force = 0.0;

      // Distance and shift function
      Real dist = sqrt(distSquared);
      info[atom1].internal.push_back(Internal(atom2, diff, dist));
      //Real rDist = 1.0 / dist;

      Real f_ij = BornSwitchValue(dist);
      Real fp_ij = BornSwitchDerivative(dist);
      Real gamma = 0.5;

      // Atom 1 variables
      int type1 = topo->atoms[atom1].type;
      Real B_i = topo->atomTypes[type1].mySCPISM->B_i;
      Real C_i = topo->atomTypes[type1].mySCPISM->C_i;
      //Real dR_vdw2 = topo->atomTypes[type1].dR_vdw2;
      Real dR_vdw2 = topo->atoms[atom1].mySCPISM->dR_vdw2;

      //**********************************
      // Part of Eq. 1
      topo->atoms[atom1].mySCPISM->sasaFrac -= B_i * f_ij * exp(-C_i * dist);
      //**********************************

      //**********************************
      // Eq. 6
      info[atom1].dRs[atom2] =
        gamma * B_i * dR_vdw2 *
        (fp_ij * exp(-C_i * dist) - C_i * f_ij * exp(-C_i * dist));
      //**********************************

      // If atom 1 is a polar H+, accumulate the derivative dR.
      if (topo->atomTypes[topo->atoms[atom1].type].mySCPISM->isHbonded == PH &&
          topo->atomTypes[topo->atoms[atom2].type].mySCPISM->isHbonded == PA &&
          (topo->atoms[atom1].residue_seq !=
           topo->atoms[atom2].residue_seq)) { // Polar
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom1].name == "HN" &&
            topo->atoms[atom2].name == "O")
          g_i = -0.378;
        else
          g_i = topo->atomTypes[type1].mySCPISM->g_i;
        Real g_j = topo->atomTypes[topo->atoms[atom2].type].mySCPISM->g_i;
        topo->atoms[atom1].mySCPISM->polarFrac += g_i * g_j * f_ij * exp(
          -E_i * dist);
        //**********************************
        // Eq. 8
        info[atom1].dRs[atom2] +=  g_i * g_j *
          (fp_ij * exp(-E_i * dist) + f_ij * exp(-E_i * dist) * -E_i);
        //**********************************
      }

      // Atom 2 variables
      int type2 = topo->atoms[atom2].type;
      B_i = topo->atomTypes[type2].mySCPISM->B_i;
      dR_vdw2 = topo->atoms[atom2].mySCPISM->dR_vdw2;
      C_i = topo->atomTypes[type2].mySCPISM->C_i;

      //**********************************
      // Part of Eq. 1
      topo->atoms[atom2].mySCPISM->sasaFrac -= B_i * f_ij * exp(-C_i * dist);
      //**********************************

      //**********************************
      // Eq. 6
      info[atom2].dRs[atom1] =  gamma * B_i * dR_vdw2 *
        (fp_ij * exp(-C_i * dist) - C_i * f_ij * exp(-C_i * dist));
      //**********************************

      // If atom 2 is a polar H+, accumulate polar
      // fraction and derivative
      if (topo->atomTypes[type2].mySCPISM->isHbonded == PH &&
          topo->atomTypes[topo->atoms[atom1].type].mySCPISM->isHbonded == PA &&
          (topo->atoms[atom2].residue_seq !=
           topo->atoms[atom1].residue_seq)) {
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom2].name == "HN" && topo->atoms[atom1].name == "O")
          g_i = -0.378;
        else g_i = topo->atomTypes[type2].mySCPISM->g_i;
        Real g_j = topo->atomTypes[type1].mySCPISM->g_i;
        topo->atoms[atom2].mySCPISM->polarFrac +=
          g_i * g_j * f_ij * exp(-E_i * dist);
        //**********************************
        // Eq. 8
        info[atom2].dRs[atom1] += g_i * g_j *
          (fp_ij * exp(-E_i * dist) + f_ij * exp(-E_i * dist) * -E_i);
        //**********************************
      }
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      // Nothing to accumulate here
    }

    static Real getEnergy(const ScalarStructure *energies) {return 0;}

    // Parsing
    static std::string getId() {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    static unsigned int getParameterSize() {return 0;}

    static CoulombBornRadiiForce make(const std::vector<Value> &);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sub Classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    class C1 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr)                    {return rr;}

      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        return cr * (1.5 - 0.5 * r * r * cr * cr);
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 1.5 * cr;}

      static Real dSmooth(Real r, Real /*c*/, Real cr) {
        return -r * cr * cr * cr;
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {
        return CoulombBornRadiiForce ::keyword;
      }

    public:
      static const std::string keyword;
    };
  public:
    class C2 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}

      static Real dKernel(Real r)                 {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr)                    {return rr;}

      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return cr * (1.875 - r *  (1.25 - 0.375 * r));
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 1.875 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * cr * cr;
        return cr * c * (1.5 * r * c - 2.5);
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) :  kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {
        return CoulombBornRadiiForce::keyword;
      }

    public:
      static const std::string keyword;
    };
  public:
    class C3 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}

      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr)                    {return rr;}

      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0625 * cr * (35.0 - r * (35.0 - r * (21.0 - 5.0 * r)));
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 2.1875 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * r * cr * cr;
        return r * cr * cr * cr * (-4.375 + c * (5.25 - 1.875 * c));
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {
        return CoulombBornRadiiForce::keyword;
      }

    public:
      static const std::string keyword;
    };
  public:
    class C4 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}

      static Real dKernel(Real r)                {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr)                    {return rr;}

      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0078125 * cr * 
          (315.0 - r * (420.0 - r * (378.0 - r * (180.0 - r * 35.0))));
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 2.4609375 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * r * cr * cr;
        return -r * cr * cr * cr *
          (6.5625 - c * (11.8125 - c * (8.4375 - c * 2.1875)));
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) :  kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {
        return CoulombBornRadiiForce::keyword;
      }

    public:
      static const std::string keyword;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    //Vector3DBlock* dR_dr;

    void resize(int newsize) {
      info.resize(newsize);
    }

    void clear() {
      info.clear();
    }

    int getSize(int i) {return info[i].internal.size();}
    void resizeDR(int i, int size) {info[i].dRs.resize(size);}
    int getInteraction(int i, int j) {return info[i].internal[j].atom2;}
    Vector3D getDiff(int i, int j) {return info[i].internal[j].diff;}
    Real getDist(int i, int j) {return info[i].internal[j].dist;}
    Real getDR(int i, int j) {return info[i].dRs[j];}

  private:
    int sw;

  public:
    Real BornSwitchValue(Real r) const {
      Real f;
      if (r > 5) return 0;
      if (sw == 0)
        return 1;
      else if (sw == 1) {
        f = (1 - .04 * r * r);
        return f * f;
      } else if (sw == 2) {
        f = (1 - 0.008 * r * r * r);
        return f * f * f;
      } else {
        //sw == 3
        f = (1 - .0016 * r * r * r * r);
        return f * f * f * f;
      }
    }

    Real BornSwitchDerivative(Real r) const {
      Real f, fr;
      if (r > 5) return 0;
      if (sw == 0)
        return 0;
      else if (sw == 1) {
        fr = -0.16 * (1 - 0.04 * r * r) * r;
        return fr;
      } else if (sw == 2) {
        fr = -9 * .008 * r * r * (1 - .008 * r * r * r) *
          (1 - .008 * r * r * r);
        return fr;
      } else {
        //sw==3
        f = (1 - .0016 * r * r * r * r);
        fr = -16 * .0016 * r * r * r * f * f * f;
        return fr;
      }
    }

    Real BornSwitchSecondDerivative(Real r) {
      Real f, frr;
      if (r > 5) return 0;
      if (sw == 0)
        return 0;
      else if (sw == 1) {
        frr = 0.0192 * r * r - 0.16;
        return frr;
      } else if (sw == 2) {
        f = (1 - 0.008 * r * r * r);
        frr = 0.003456 * f * r * r * r * r - 0.144 * f * f * r;
        return frr;
      } else if (sw == 3) {
        //sw==3 function
        f = (1 - .0016 * r * r * r * r);
        frr = 0.00049152 * f * f * r * r * r * r * r * r - 0.0768 * f * f * f *
              r * r;
        return frr;
      } else {
        f = (1 - 0.04 * r * r);
        frr = 0.0384 * f * r * r - 0.24 * f * f;
        return frr;
      }
    }
  };

  //____ INLINES
}
#endif /* COULOMBFORCE_H */
