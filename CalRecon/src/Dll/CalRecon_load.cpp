/** 
* @file CalRecon_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
* $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

DECLARE_FACTORY_ENTRIES(CalRecon) {
    DECLARE_ALGORITHM( CalXtalRecAlg);
    DECLARE_ALGORITHM( CalClustersAlg);
    DECLARE_ALGORITHM( CalDisplay);
    DLL_DECL_TOOL( SingleClusterTool );
    DLL_DECL_TOOL( LastLayerCorrTool );
    DLL_DECL_TOOL( ProfileTool );
} 




