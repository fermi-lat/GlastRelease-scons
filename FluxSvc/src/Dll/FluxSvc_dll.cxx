//====================================================================
//  FluxSvc_dll.cpp
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
#include "GaudiKernel/DllMain.icpp"

void GaudiDll::initialize(void* /* hinstDLL */ )    {
}

void GaudiDll::finalize(void* /* hinstDLL */ )      {
}

extern void FluxSvc_load();

#include "GaudiKernel/FactoryTable.h"

extern "C" FactoryTable::EntryList* getFactoryEntries() {
    static bool first = true;
    if ( first ) {
        FluxSvc_load();
        first = false;
    }
    return FactoryTable::instance()->getEntries();
} 


