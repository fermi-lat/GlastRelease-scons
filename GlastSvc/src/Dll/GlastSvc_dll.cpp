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
//  Author     : H. Gillespie
//  Created    : 1 Aug 00
//  Changes    : 
//               Version copied from Glast, and addopted for GlastSvc
//
//====================================================================

// DllMain entry point
#include "Gaudi/System/DllMain.icpp"
#include <iostream>
void GaudiDll::initialize(void*) 
{
}

void GaudiDll::finalize(void*) 
{
}
extern void GlastSvc_load();
#include "Gaudi/Kernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {  // Don't call for UNIX
    GlastSvc_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 



