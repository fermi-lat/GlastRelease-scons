//====================================================================
//  GuiSvc_load.cpp
//--------------------------------------------------------------------
//
//  Package    : GuiSvc
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

void GuiSvc_load() {
  DLL_DECL_SERVICE( GuiSvc );

  // declare the algorithms needed to implement the event loop

}

extern "C" void GuiSvc_loadRef() {
  GuiSvc_load();
}
