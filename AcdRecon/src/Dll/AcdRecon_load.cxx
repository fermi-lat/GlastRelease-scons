/** 
* @file AcdRecon_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "../AcdReconAlg.h"
// #include "../AcdDisplay.h"   // need to separated out AcdDisplay.h 

DECLARE_ALGORITHM_FACTORY( AcdReconAlg );
// DECLARE_ALGORITHM_FACTORY( AcdDisplay );  // temp done in AcdDisplay.cxx

DECLARE_FACTORY_ENTRIES( AcdRecon ) {
    DECLARE_ALGORITHM( AcdReconAlg );
    DECLARE_ALGORITHM( AcdDisplay );
} 
