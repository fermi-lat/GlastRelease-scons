/** 
 * @file ntupleWriterSvc_dll.cxx
 *
 * @brief Implementation of the DllMain routine.  
 * The DLL initialization must be done seperately for each DLL. 
 * 
 * @author H. Gillespie
 * 
 * $Header$
 */

// DllMain entry point
#include "GaudiKernel/DllMain.icpp"
#include <iostream>
void GaudiDll::initialize(void*) 
{
}

void GaudiDll::finalize(void*) 
{
}
extern void ntupleWriterSvc_load();
#include "GaudiKernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {  // Don't call for UNIX
    ntupleWriterSvc_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 

void FATAL(const char * msg) {
    std::cerr << "FATAL error from ntupleWriterSvc DLL: " << msg << std::endl;
}

void WARNING(const char * t){ std::cerr << "WARNING from ntupleWriterSvc DLL: " << t << std::endl; }
