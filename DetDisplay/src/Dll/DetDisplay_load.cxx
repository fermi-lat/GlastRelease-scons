/** 
* @file DetDisplay_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/
#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_FACTORY_ENTRIES(DetDisplay) {
    DECLARE_TOOL( StripDisplay );
    DECLARE_TOOL( MCdisplay );
    DECLARE_TOOL( DetectorDisplay );
} 
