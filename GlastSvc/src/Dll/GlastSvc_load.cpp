/** 
* @file GlastSvc_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(GlastSvc) {
    DECLARE_SERVICE( GlastDetSvc);
    DECLARE_SERVICE( GlastEventSelector );
    DECLARE_SERVICE( EventCnvSvc );
    
    DECLARE_CONVERTER( EventCnv );
    DECLARE_CONVERTER( MCEventCnv );
   //THB DECLARE_CONVERTER( DigiEventCnv );
} 
