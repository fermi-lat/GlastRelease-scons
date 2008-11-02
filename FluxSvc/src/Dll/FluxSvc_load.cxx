/** 
* @file FluxSvc_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/
#include "GaudiKernel/DeclareFactoryEntries.h"


DECLARE_FACTORY_ENTRIES(FluxSvc) {
    DECLARE_SERVICE( FluxSvc );
    DECLARE_SERVICE( OrbitSvc );
    DECLARE_ALGORITHM( FluxAlg );
    DECLARE_ALGORITHM( ExposureAlg );
    DECLARE_ALGORITHM( PointInfoAlg );
    DECLARE_TOOL( FluxSvcRandom );
}

