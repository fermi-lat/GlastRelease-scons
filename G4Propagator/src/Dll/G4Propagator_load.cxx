/** 
* @file G4Propagator_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(G4Propagator) {
    DECLARE_TOOL( G4PropagatorTool );
    DECLARE_TOOL( G4PropagationTool );
} 
