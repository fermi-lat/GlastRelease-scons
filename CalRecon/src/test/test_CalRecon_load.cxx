//====================================================================
//
//  Description: Implementation of <Package>_load routine.
//               This routine is needed for forcing the linker
//               to load all the components of the library. 
//
//====================================================================

#include "GaudiKernel/DeclareFactoryEntries.h"

//! Load all  services: 
void test_CalRecon_load() {
    DECLARE_ALGORITHM( test_CalRecon );
} 

extern "C" void test_CalRecon_loadRef()    {
    test_CalRecon_load();
}
