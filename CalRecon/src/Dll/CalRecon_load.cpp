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

//! Load all  services: 
void CalRecon_load() {
    
//    DLL_DECL_SERVICE( CalGeometrySvc );
    
//    DLL_DECL_ALGORITHM( CalRecLogsAlg );
//    DLL_DECL_ALGORITHM( CalDigiLogsAlg );
//     DLL_DECL_ALGORITHM( CalIRFAlg );
    DLL_DECL_ALGORITHM( CalXtalRecAlg );
    DLL_DECL_ALGORITHM( CalClustersAlg );
//    DLL_DECL_ALGORITHM( CalNtupleAlg );
//    DLL_DECL_ALGORITHM( CalDisplay );
} 

extern "C" void CalRecon_loadRef()    {
    CalRecon_load();
}

