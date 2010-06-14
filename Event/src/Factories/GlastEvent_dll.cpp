//====================================================================
//  GlastEventFactories_dll.cpp
//--------------------------------------------------------------------
//
//  Package    : Gaudi
//
//  Description: Implementation of DllMain routine.
//               The DLL initialisation must be done seperately for 
//               each DLL. 
//
//  Author     : M.Frank
//  Created    : 13/1/99
//  Changes    : 
//
//====================================================================

// DllMain entry point
#include "Gaudi/System/DllMain.icpp"
#include "Gaudi/Kernel/FactoryTable.h"
#include <iostream>
extern void GlastEventFactories_load();

void GaudiDll::initialize(void* /*hinstDLL*/)    {
}

void GaudiDll::finalize(void* /*hinstDLL*/)    {
}

extern "C" FactoryTable::EntryList* getFactoryEntries() {
  static bool first = true;
  if ( first )    {
    GlastEventFactories_load();
    first = false;
  }
  return FactoryTable::instance()->getEntries();
} 


void    FATAL ( const char* c ) 
{
    std::cerr << "Glast based FATAL error: " << c << std::endl;
}

void    WARNING ( const char* c ) 
{
    std::cerr << "Glast based WARNING error: " << c << std::endl;
}

