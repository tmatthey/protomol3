#ifndef MODIFIEROUTPUTINTERMEDIATEPOINTS_H
#define MODIFIEROUTPUTINTERMEDIATEPOINTS_H

#include <protomol/modifier/Modifier.h>
#include <src/force/bonded/HarmNMRestSystemForce.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <protomol/io/PDBWriter.h>

namespace ProtoMol {

   class Integrator;

   //  _________________________________________________________ ModifierOutputIntermediatePoints
   
   class ModifierOutputIntermediatePoints : public Modifier {

       typedef HarmNMRestSystemForce<VacuumBoundaryConditions> HarmNMForce;

        //  ----------------------------------------------------------------  //
        //  Constructors, destructors, assignment
        //  ----------------------------------------------------------------  //
 
       public:
          ModifierOutputIntermediatePoints(const Integrator *i, int n_p);
          
          ~ModifierOutputIntermediatePoints();

        //  ----------------------------------------------------------------  //
        //  From class Modifier
        //  ----------------------------------------------------------------  //

       private:
          virtual void doExecute(Integrator* i);
          virtual void doInitialize();

       public:
          virtual bool isInternal() const { return( false ); }
          virtual std::string getIdNoAlias() const {return "OutputIntermediatePoints";}


 
        //  ------------------------------------------------------------------ //
        //  New function of the class ModifierOutputIntermediatePoints         //
        //  ------------------------------------------------------------------ //
 
       protected:
          virtual bool GetForcePointer(std::string force_name);
          virtual bool IdentifyForce(string force_name, const Integrator *integrator);
          string NextIntermediatePointFilename();
          Real getTimestep() const;

        // ------------------------------------------------------------------- //
        // class members/variables                                             //
        // ------------------------------------------------------------------- //
       private:
          const Integrator *myIntegrator;
 
          HarmNMForce *fptr;

          //Need a PDBWriter
          PDBWriter pdbWriter;

          int next_point;
          std::string pts_filename;

   };

}

#endif
