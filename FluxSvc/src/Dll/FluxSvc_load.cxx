//====================================================================
//  FluxSvc_load.cpp
//--------------------------------------------------------------------
//
//  Package    : FluxSvc
//
//  Description: Implementation of <Package>_load routine. This routine 
//               is needed for forcing the linker to load all the components
//               of the library.. 
//
//====================================================================
#include "GaudiKernel/ICnvFactory.h"
#include "GaudiKernel/ISvcFactory.h"
#include "GaudiKernel/IAlgFactory.h"

#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();

void FluxSvc_load() {
  DLL_DECL_SERVICE( FluxSvc );


}

extern "C" void FluxSvc_loadRef() {
  FluxSvc_load();
}
