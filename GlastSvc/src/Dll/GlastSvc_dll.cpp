//====================================================================
//  GlastSvc_dll.cpp
//--------------------------------------------------------------------
//
//  Package    : GlastSvc
//
//  Description: Implementation of DllMain routine.
//               The DLL initialisation must be done seperately for 
//               each DLL. 
//
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
extern void GlastSvc_load();
#include "GaudiKernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {  // Don't call for UNIX
    GlastSvc_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 


