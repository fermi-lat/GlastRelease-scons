//====================================================================
//  GlastSvc_load.cpp
//--------------------------------------------------------------------
//
//  Package    : Gaudi/System
//
//  Description: Implementation of <Package>_load routine.
//               This routine is needed for forcing the linker
//               to load all the components of the library. 
//
//====================================================================

#include "GaudiKernel/DeclareFactoryEntries.h"


DECLARE_FACTORY_ENTRIES(McToolBox) 
{
    DECLARE_ALGORITHM( McToolBoxAlg            );

    DECLARE_TOOL(      McBuildRelTablesTool    );
    DECLARE_TOOL(      McGetEventInfoTool      );
    DECLARE_TOOL(      McGetTrackInfoTool      );
} 

