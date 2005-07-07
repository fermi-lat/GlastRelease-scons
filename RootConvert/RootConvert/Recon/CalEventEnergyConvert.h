#ifndef RootConvert_CalEventEnergyConvert_H
#define RootConvert_CalEventEnergyConvert_H 1

#include <Event/Recon/CalRecon/CalEventEnergy.h>   
#include <reconRootData/CalEventEnergy.h>

namespace RootPersistence {

void convert( const Event::CalCorToolResult &, CalCorToolResult & ) ;
void convert( const CalCorToolResult &, Event::CalCorToolResult & ) ;
void convert( const Event::CalEventEnergy &, CalEventEnergy & ) ;
void convert( const CalEventEnergy &, Event::CalEventEnergy & ) ;

}

#endif












