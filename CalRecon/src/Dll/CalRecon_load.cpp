/** 
* @file CalRecon_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
* $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

DECLARE_FACTORY_ENTRIES(CalRecon) {
    DECLARE_ALGORITHM( CalXtalRecAlg);
    DECLARE_ALGORITHM( CalClustersAlg);
    DECLARE_ALGORITHM( CalDisplay);
    DECLARE_TOOL( SingleClusterTool );
    DECLARE_TOOL( FuzzyClusterTool );
    DECLARE_TOOL( SimpleClusterTool );
    DECLARE_TOOL( LastLayerCorrTool );
    DECLARE_TOOL( ProfileTool );
    DECLARE_TOOL( CalValsCorrTool );
} 




