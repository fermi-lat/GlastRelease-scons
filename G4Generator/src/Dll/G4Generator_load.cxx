/** 
* @file G4Generator_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(G4Generator) {
    DECLARE_ALGORITHM( G4Generator);
} 
