/*! \file GlastEventFactories_dll.cpp
\brief based upon LHCbEventFactories_dll.cpp by M. Frank available within the LHCbEvent package

  Implementation of DllMain routine.
  The DLL initialisation must be done seperately for each DLL.
*/

// DllMain entry point
#include "GaudiKernel/DllMain.icpp"
#include "GaudiKernel/FactoryTable.h"
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

