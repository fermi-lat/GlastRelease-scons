/** 
* @file TkrUtil_load.cxx
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(TkrUtil) {

    DECLARE_ALGORITHM ( TkrFillTDInfoAlg  );
    
    DECLARE_SERVICE( TkrFailureModeSvc   );
    DECLARE_SERVICE( TkrBadStripsSvc     );
    DECLARE_SERVICE( TkrGeometrySvc      );
    DECLARE_SERVICE( TkrAlignmentSvc     );
    DECLARE_SERVICE( TkrSplitsSvc        );
    DECLARE_SERVICE( TkrToTSvc           );

    DECLARE_TOOL(    TkrQueryClustersTool);
    DECLARE_TOOL(    TkrMakeClustersTool );
    DECLARE_TOOL(    TkrEnergyTool       );

    DECLARE_ALGORITHM ( TkrCalibAlg      );

    DECLARE_TOOL(    TkrFlagHitsTool     );
    DECLARE_TOOL(    TkrGhostTool        );
    DECLARE_TOOL(    TkrDiagnosticTool   );
    DECLARE_TOOL(    TkrMapTool          );
    DECLARE_TOOL(    TkrTrackVecTool     );
    DECLARE_TOOL(    TkrReasonsTool      );
    DECLARE_TOOL(    TkrHitTruncationTool);
} 



