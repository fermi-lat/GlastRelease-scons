/** 
* @file AcdDigi_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(AcdDigi) {
    DECLARE_ALGORITHM( AcdDigiAlg );
    DECLARE_ALGORITHM( AcdDigiMcIntHitAlg );
    DECLARE_ALGORITHM( AcdDigiOrgAlg );
} 

