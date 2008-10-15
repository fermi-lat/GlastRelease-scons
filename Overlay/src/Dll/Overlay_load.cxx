/** 
* @file RootIo_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(Overlay) {
    DECLARE_ALGORITHM( overlayRootReaderAlg );
    DECLARE_TOOL( BackgroundSelectTool );
    DECLARE_TOOL( McIlwain_L_Tool );
}
  
