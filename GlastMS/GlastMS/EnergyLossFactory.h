/** @file EnergyLossFactory.h
 **
 **
 **  $Header$
 **
 */
#ifndef EnergyLossFactory_h
#define EnergyLossFactory_h

// Interface class definition for energy loss
// For Geant4 this is a both a "continuous" and "discreete" process
#include "G4VContinuousDiscreteProcess.hh"

namespace GlastMS {
    /** @class EnergyLossFactory
        @brief simple class to create a G4 process for energy loss
               including both ionization and bremstrahlung 
               (and pair production for muons!)
    **/
    class EnergyLossFactory {
    public:
        enum ELossVersion  {CURRENT,    RELEASE52};
        enum ELossType     {IONIZATION, BREMSSTRAHLUNG, PAIRPRODUCTION};
        enum Particle      {ELECTRON,   POSITRON,   MUON};

        /** contstructor chooses which one
        @param release[CURRENT]   to get it from the current release of G4
               release[RELEASE52] to get it from G4 version 5.2
        */
        EnergyLossFactory(ELossVersion release=CURRENT):m_release(release){}

        /**run the factory: this is a function operator
        @param part[ELECTRON,POSITRON] for electrons and positrons
               part[MUON]              for muons
        @param type[IONIZATION]        ionization energy loss
               type[BREMSSTRAHLUNG]    Bremsstrahlung losses (which result in delta rays)
               type[PAIRPRODUCTION]    For muons only
        */
     G4VContinuousDiscreteProcess* operator()(Particle part = ELECTRON, ELossType type = IONIZATION);

    private:
       ELossVersion m_release;
    };

}// namespace

#endif
