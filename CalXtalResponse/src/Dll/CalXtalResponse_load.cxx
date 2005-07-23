/** 
 * @file CalXtalResponse_load.cpp
 * @brief This is needed for forcing the linker to load all components
 * of the library.
 *
 *  $Header$
 */

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

DECLARE_FACTORY_ENTRIES(CalXtalResponse) {
  
  DECLARE_SERVICE( CalCalibSvc );
  
  DECLARE_ALGORITHM( CalXtalRecAlg);

  DLL_DECL_TOOL( XtalDigiTool );
  DLL_DECL_TOOL( XtalRecTool );
  
  DLL_DECL_TOOL( CalXtalRespRandom );
} 
