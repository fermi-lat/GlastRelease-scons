//====================================================================
//  GuiSvc_dll.cpp
//--------------------------------------------------------------------
//
//  Package    : Gaudi/System
//
//  Description: Implementation of DllMain routine.
//               The DLL initialisation must be done seperately for 
//               each DLL. 
//
//  Author     : Toby Burnett
//
//====================================================================

// DllMain entry point
#include "Gaudi/System/DllMain.icpp"

void GaudiDll::initialize(void* /* hinstDLL */ )    {
}

void GaudiDll::finalize(void* /* hinstDLL */ )      {
}

extern void GuiSvc_load();
#include "Gaudi/Kernel/FactoryTable.cpp"
extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first ) {
    GuiSvc_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 


