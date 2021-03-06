/** 
* @file EbfWriter_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

// There should be one entry for each component included in 
// the library for this package.
DECLARE_FACTORY_ENTRIES(EbfWriter) {
    DECLARE_ALGORITHM(EbfWriter);
      DECLARE_ALGORITHM(FesExposure);
} 
