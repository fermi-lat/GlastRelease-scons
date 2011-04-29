/** @file EnergyLossFactory.cxx

   $Header$
   


*/
#include "GlastMS/EnergyLossFactory.h"
#include "G4ExceptionHandler.hh"
#include "G4eIonisation.hh"
//#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung.hh"
//#include "G4eBremsstrahlung52.hh"
#include "G4MuIonisation.hh"
//#include "G4MuIonisation52.hh"
#include "G4MuBremsstrahlung.hh"
//#include "G4MuBremsstrahlung52.hh"
#include "G4MuPairProduction.hh"
//#include "G4MuPairProduction52.hh"

G4VContinuousDiscreteProcess* GlastMS::EnergyLossFactory::operator()(Particle part, ELossType type)
{
    // Pointer to what we are trying to make...
    G4VContinuousDiscreteProcess* process = 0;

    // First determine from which version of code we are going to look at 
    if( m_release == CURRENT)
    {
        // What kind of particle are we dealing with?
        if (part == ELECTRON || part == POSITRON)
        {
            // Now look at what type of Energy Loss we are creating. 
            if   (type == IONIZATION) process = new G4eIonisation();
            else                      process = new G4eBremsstrahlung();
        }
        else  // assume Muon
        {
            // Now look at what type of Energy Loss we are creating. 
            if      (type == IONIZATION)     process = new G4MuIonisation();
            else if (type == BREMSSTRAHLUNG) process = new G4MuBremsstrahlung();
            else                             process = new G4MuPairProduction();
        }
    }
    else if (m_release == RELEASE52) // NOT AVAILABLE ANYMORE
    {
        if (part == ELECTRON || part == POSITRON)
        {
            // Now look at what type of Energy Loss we are creating. 
            if   (type == IONIZATION) process = new G4eIonisation();
            else                      process = new G4eBremsstrahlung();
        }
        else  // assume Muon
        {
            // Now look at what type of Energy Loss we are creating. 
            if      (type == IONIZATION)     process = new G4MuIonisation();
            else if (type == BREMSSTRAHLUNG) process = new G4MuBremsstrahlung();
            else                             process = new G4MuPairProduction();
        }
    }
    else
    {
        // Best not to get to here...
        G4Exception("Illegal Releast type in EnergyLossFactory!");
    }

    return process;
}
