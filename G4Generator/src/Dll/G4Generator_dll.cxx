//====================================================================
//
//====================================================================

// DllMain entry point
#include "GaudiKernel/DllMain.icpp"
#include <iostream>
void GaudiDll::initialize(void*) 
{
}

void GaudiDll::finalize(void*) 
{
}
extern void G4Generator_load();
#include "GaudiKernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {  // Don't call for UNIX
    G4Generator_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 


