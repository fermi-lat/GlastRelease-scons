/** 
* @file RootIo_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(Overlay) {
    DECLARE_ALGORITHM( AcdOverlayMergeAlg );
    DECLARE_ALGORITHM( CalOverlayMergeAlg );
    DECLARE_ALGORITHM( TkrOverlayMergeAlg );
    DECLARE_ALGORITHM( GemOverlayMergeAlg );
    DECLARE_ALGORITHM( DigiToOverlayAlg );
    DECLARE_SERVICE( OverlayInputSvc );
    DECLARE_SERVICE( OverlayOutputSvc );
    DECLARE_TOOL( McIlwain_L_Tool );
    DECLARE_TOOL( EventToOverlayTool );
    DECLARE_TOOL( CalXtalToOverlayTool );
    DECLARE_TOOL( TkrDigiToOverlayTool );
    DECLARE_TOOL( AcdHitToOverlayTool );
    DECLARE_TOOL( GemToOverlayTool );
    DECLARE_TOOL( OverlayRandom);
    DECLARE_CONVERTER( SrcOverlayCnv );
    DECLARE_CONVERTER( EventOverlayCnv );
    DECLARE_CONVERTER( TkrOverlayCnv );
    DECLARE_CONVERTER( CalOverlayCnv );
    DECLARE_CONVERTER( AcdOverlayCnv );
    DECLARE_CONVERTER( GemOverlayCnv );
//    DECLARE_CONVERTER( LdfGemCnv );
}
  
