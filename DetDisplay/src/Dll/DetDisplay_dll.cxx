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
extern void DetDisplay_load();
#include "GaudiKernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {  // Don't call for UNIX
    DetDisplay_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 

void FATAL(const char * msg) {
    std::cerr << "Fatal error from DetDisplay DLL: " << msg << std::endl;
}

void WARNING(const char * msg) {
    std::cerr << "Warning message from DetDisplay DLL:" << msg << std::endl;
}

