/** 
* @file CalRecon_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES( CalRecon ) 
{
    DECLARE_SERVICE(   CalReconSvc             );
    DECLARE_ALGORITHM( CalClustersAlg          );
    DECLARE_ALGORITHM( CalMipFinderAlg         );
    DECLARE_ALGORITHM( CalEventEnergyAlg       );
    DECLARE_ALGORITHM( CalDisplay              );
    DECLARE_ALGORITHM( PropertiesCheckAlg      );
    DECLARE_TOOL(      CalSingleClusteringTool );
    DECLARE_TOOL(      CalSimpleClusteringTool );
    DECLARE_TOOL(      StdMipFindingTool       );
    DECLARE_TOOL(      CalRawEnergyTool        );
    DECLARE_TOOL(      CalTkrLikelihoodTool    );
    DECLARE_TOOL(      CalLastLayerLikelihoodTool );
    DECLARE_TOOL(      CalFullProfileTool      );
    DECLARE_TOOL(      CalProfileTool          );
    DECLARE_TOOL(      CalValsCorrTool         );
    DECLARE_TOOL(      CalTransvOffsetTool     );
} 




