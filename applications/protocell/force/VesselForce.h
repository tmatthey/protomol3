/* -*- c++ -*- */
#ifndef VESSELFORCE_H
#define VESSELFORCE_H

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>

#include <string>

#include <protomol/base/Report.h>
using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ VesselForce

  template<class TBoundaryConditions>
  class VesselForce : public SystemForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    VesselForce() : mySigma(0.0), myEpsilon(0.0), 
                    myDiameter(0.0), myGranularity(0.0) {}
    VesselForce( double sig, double eps, double dia, double gran ) : mySigma(sig),
        myEpsilon(eps), myDiameter(dia), myGranularity(gran) {}
    virtual ~VesselForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class VesselForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void calcInteraction(const GenericTopology *topo,
                  const TBoundaryConditions &boundary, const int atomIndex,
                  const Vector3DBlock *positions, Vector3DBlock *forces,
                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions, Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const {return "Vessel";}

    virtual unsigned int numberOfBlocks(const GenericTopology *topo,
                                        const Vector3DBlock *pos);

  private:
    virtual Force *doMake(const std::vector<Value> &) const;// {
    //  return new VesselForce();
    //}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return getKeyword();}
    virtual void getParameters(std::vector<Parameter> &) const;// {}
    virtual unsigned int getParameterSize() const;//{return 2;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    double mySigma, myEpsilon, myDiameter, myGranularity;

  };

  //____ INLINES

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::evaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      ((SemiGenericTopology<TBoundaryConditions> &)(*topo)).
        boundaryConditions;

    const unsigned int count = positions->size();

    for (unsigned int i = 0; i < count; i++){
        calcInteraction(topo, boundary, i, positions, forces, energies);
    }
    
  }


  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::parallelEvaluate(
    const GenericTopology *topo, const Vector3DBlock *positions,
    Vector3DBlock *forces, ScalarStructure *energies) {
    const TBoundaryConditions &boundary =
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions> &>(*topo)).
        boundaryConditions;

    unsigned int n = positions->size();
    unsigned int count = numberOfBlocks(topo, positions);

    for (unsigned int i = 0; i < count; i++)
      if (Parallel::next()) {
        int to = (n * (i + 1)) / count;
        if (to > static_cast<int>(n))
          to = n;
        int from = (n * i) / count;
        for (int j = from; j < to; j++){
            calcInteraction(topo, boundary, i, positions, forces, energies);
        }
      }

  }

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::calcInteraction(
    const GenericTopology *topo,
    const TBoundaryConditions &boundary, const int atomIndex,
    const Vector3DBlock *positions, Vector3DBlock *forces,
    ScalarStructure *energies) {

    //get atom type
    const int type = topo->atoms[atomIndex].type;
    //Get parameters
    Real sigma = topo->atomTypes[type].sigma;
    Real epsilon = topo->atomTypes[type].epsilon;
    //Real vesselSigma = 1.15;
    //Real vesselEpsilon = -0.01;//-0.000511;
    Real r_ij = sigma + mySigma;
    Real e_ij = sqrt(epsilon * myEpsilon);
    Real A = power<12>(r_ij) * e_ij;
    Real B = 2 * power<6>(r_ij) * e_ij;

    //report << hint << "i " << atomIndex << " name " << topo->atoms[atomIndex].name <<  " A " << A << " B " << B << " sigma " << sigma << " epsilon "<< epsilon << endr;
    
    //Calculate SCE-Vessel force
    // Distance squared, diff and offs vector
    Real distSquared(0.0);

    Vector3D diff(0.0,0.0,0.0);

    Vector3D offs(0.0,0.0,0.0);

    //get MINIMAL difference vector, removes PBC extents
    Vector3D diff1 = boundary.minimalDifference((*positions)[atomIndex], offs );

    //assume cylinder diameter
    const Real cylinderDia = myDiameter;

    // find nearest 'side' of cylinder (either add or subtract the diameter)

    //Y- component first
    //ydiff will be the distance to cylinder, yOffs the cylinder position
    Real ydiff, yOffs;

    //smallest
    if(fabs(diff1[1] - cylinderDia) < fabs(diff1[1] + cylinderDia)){
      ydiff = diff1[1] - cylinderDia;
      yOffs = -cylinderDia;
    }else{
      ydiff = diff1[1] + cylinderDia;
      yOffs = cylinderDia;
    }

    //then Z-component
    Real zdiff, zOffs;

    if(fabs(diff1[2] - cylinderDia) < fabs(diff1[2] + cylinderDia)){
      zdiff = diff1[2] - cylinderDia;
      zOffs = -cylinderDia;
    }else{
      zdiff = diff1[2] + cylinderDia;
      zOffs = cylinderDia;
    }

    //find closest, Y or Z and change the diff vector (direction of force) and distSquared to represent it
    //diff zero here
    //
    if(fabs(ydiff) < fabs(zdiff)){ //added fabs 1/2/10, also changed to ydiff, zdiff below
      diff[1] = ydiff;//yOffs;
      distSquared = ydiff * ydiff;
    }else{
      diff[2] = zdiff;//zOffs;
      distSquared = zdiff * zdiff;
    }

    //figure out granularity of vessel, to represent cells
    const double vessel_granularity = myGranularity;

    //find last "cell" and next cell
    double last_vessel_cell = diff1[0] - floor( diff1[0] / vessel_granularity ) * vessel_granularity;

    //double next_vessel_cell = vessel_granularity - last_vessel_cell;

    //find smallest
    //if( fabs(last_vessel_cell) < fabs(next_vessel_cell)){
      //if(diff1[0] > 0.0 )
        diff[0] = last_vessel_cell; //orig -
      //else
      //  diff[0] = last_vessel_cell;

    //}else{
      //if(diff1[0] > 0.0 )
        //diff[0] = next_vessel_cell; //orig +
      //else
      //  diff[0] = -next_vessel_cell;
    //}

    //get distance to vessel cell
    distSquared = diff.normSquared();

    //report << hint << "Distance " << distSquared << endr;
    
    //Test!
    if( distSquared ){
        Real rDistSquared = 1.0 / distSquared;

        // Calculate force on atom due to vessel.
        Real r6 = rDistSquared * rDistSquared * rDistSquared;
        Real r12 = r6 * r6;
        Real r6B = B * r6;
        Real r12A = A * r12;
        Real energy = r12A - r6B;
        Real force = 12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared;

        //report << hint << "Force " << force << endr;

        //add it
        Vector3D fiv(diff * force);

        // Add to the total force.
        (*forces)[atomIndex] -= fiv;

        // Add energy
        (*energies)[ScalarStructure::LENNARDJONES] += energy; //or OTHER?
    }
    
    // Add virial
    /*if (energies->virial())
      energies->addVirial(force1, r12);*/
  }

  template<class TBoundaryConditions>
   unsigned int VesselForce<TBoundaryConditions>::numberOfBlocks(
    const GenericTopology *topo, const Vector3DBlock *positions) {
    return std::min(Parallel::getAvailableNum(),
                    static_cast<int>(positions->size()));
  }

  template<class TBoundaryConditions>
   inline Force* VesselForce<TBoundaryConditions>::doMake(const std::vector<Value> & values) const {
        return new VesselForce(values[0], values[1], values[2], values[3]);

  }

  template<class TBoundaryConditions>
  inline void VesselForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const {
      parameters.push_back(Parameter("-sigma", Value(mySigma)));
      //Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
      //         80, Text("Bulk solvent dielectric"))
      parameters.push_back(Parameter("-epsilon", Value(myEpsilon)));
      parameters.push_back(Parameter("-diameter", Value(myDiameter)));
      parameters.push_back(Parameter("-granularity", Value(myGranularity)));
  }

  template<class TBoundaryConditions>
  inline unsigned int VesselForce<TBoundaryConditions>::getParameterSize() const {
      return 4;
  }

}

#endif /* VESSELFORCE_H */
