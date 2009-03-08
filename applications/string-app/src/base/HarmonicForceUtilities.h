/*  -*- c++ -*-  */
#ifndef HARMONICFORCEUTILITIES_H
#define HARMONICFORCEUTILITIES_H

#include <src/force/bonded/HarmNMRestSystemForce.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/integrator/Integrator.h>
#include <string>


namespace ProtoMol {

  typedef HarmNMRestSystemForce<VacuumBoundaryConditions> HarmNMForce;    

  HarmNMForce *GetHarmNMForcePointer(std::string force_name, const Integrator *integrator);

  HarmNMForce *IdentifyHarmNMForce(std::string force_name, const Integrator *integrator);

}

#endif /* HARMONICFORCEUTILITIES_H */
