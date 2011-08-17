/*
* @file GlastClassify_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(GlastClassify) {
    DECLARE_ALGORITHM(ClassifyAlg);
    DECLARE_TOOL(ClassifyTool);
} 



