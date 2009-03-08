#ifndef MODIFIERHARMNMRESTFORCE_H
#define MODIFIERHARMNMRESTFORCE_H

#include <protomol/modifier/Modifier.h>
#include <src/force/bonded/HarmNMRestSystemForce.h>
#include <protomol/topology/VacuumBoundaryConditions.h>


namespace ProtoMol {

   class Integrator;

   //  _________________________________________________________ ModifierHarmNMRestForce

   class ModifierHarmNMRestForce : public Modifier {
       typedef HarmNMRestSystemForce<VacuumBoundaryConditions> HarmNMForce;

        //  ----------------------------------------------------------------  //
        //  Constructors, destructors, assignment
        //  ----------------------------------------------------------------  //

       public:
          ModifierHarmNMRestForce(const Integrator *i);

        //  ----------------------------------------------------------------  //
        //  From class Modifier
        //  ----------------------------------------------------------------  //

       private:
          virtual void doExecute(Integrator* i);
          virtual void doInitialize();

       public:
          virtual bool isInternal() const { return( false ); }
          virtual std::string getIdNoAlias() const {return "HarmNM";}

        //  ------------------------------------------------------------------ //
        //  New function of the class ModifierHarmNMRestForce                  //
        //  ------------------------------------------------------------------ // 
        
       protected:
          virtual Real getTimestep() const;     
          virtual bool GetForcePointer(std::string force_name);
          virtual bool IdentifyForce(string force_name, const Integrator *integrator);
            

        // ------------------------------------------------------------------- //
        // class members/variables                                             //
        // ------------------------------------------------------------------- //
       private:
          const Integrator *myIntegrator;

          HarmNMForce *fptr;


   };
}

#endif /* MODIFIERHARMNMRESTFORCE_H */
