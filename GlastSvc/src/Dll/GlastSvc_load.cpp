// $Header$
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

#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"


#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_OBJECT(x)     extern const IFactory& x##Factory; x##Factory.addRef();

//! Load all GlastSvc services: 
void GlastSvc_load() {
    // GlastSvc services
    DLL_DECL_SERVICE( GlastDetSvc );
    DLL_DECL_SERVICE( GlastEventSelector );
    DLL_DECL_SERVICE( EventCnvSvc );
    //DLL_DECL_CONVERTER( MCACDHitCnv );
    DLL_DECL_CONVERTER( EventCnv );
    DLL_DECL_CONVERTER( MCEventCnv );
    DLL_DECL_CONVERTER( DigiEventCnv );

    // DLL_DECL_CONVERTER( MCSiLayerCnv );
    //DLL_DECL_CONVERTER( MCCalorimeterHitCnv );

    DLL_DECL_CONVERTER( LdGlastDataCnv );

} 

extern "C" void GlastSvc_loadRef()    {
  GlastSvc_load();
}

