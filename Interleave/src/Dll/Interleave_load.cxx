/** 
* @file Interleave_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(Interleave) {
    DECLARE_TOOL( RegisterSampledBackground );
    DECLARE_TOOL( BkgndTupleSelectTool );
    DECLARE_ALGORITHM( InterleaveAlg);
    DECLARE_ALGORITHM( FilterFailAlg);
} 

