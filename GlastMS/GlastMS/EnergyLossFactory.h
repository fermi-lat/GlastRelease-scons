/** @file EnergyLossFactory.h
 **
 **
 **  $Header$
 **
 */
#ifndef EnergyLossFactory_h
#define EnergyLossFactory_h

// Interface class definition for energy loss
#include "G4VContinuousDiscreteProcess.hh"

namespace GlastMS {
    /** @class EnergyLossFactory
        @brief simple class to create a G4 process for energy loss
               including both ionization and bremstrahlung
    **/
    class EnergyLossFactory {
    public:
        enum ELossVersion  {CURRENT,    RELEASE52};
        enum ELossType     {IONIZATION, BREMSTRAHLUNG};
        enum Particle      {ELECTRON,   POSITRON,   MUON};

        /** contstructor chooses which one
        @param type[NATIVE] get it from the installed G4 
        */
        EnergyLossFactory(ELossVersion release=CURRENT):m_release(release){}

        /**run the factory: this is a function operator
        */
     G4VContinuousDiscreteProcess* operator()(Particle part = ELECTRON, ELossType type = IONIZATION);

    private:
       ELossVersion m_release;
    };

}// namespace

#endif