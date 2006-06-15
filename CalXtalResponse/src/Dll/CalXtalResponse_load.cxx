//  $Header$
/** 
 * @file
 * @author Zach Fewtrell
 * @brief This is needed for forcing the linker to load all components
 * of the library.
 *
 */

#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/IToolFactory.h"

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

#define DLL_DECL_TOOL(x)       extern const IToolFactory& x##Factory; x##Factory.addRef();

DECLARE_FACTORY_ENTRIES(CalXtalResponse) {
  
  DECLARE_SERVICE( CalCalibSvc );
  DECLARE_SERVICE( CalFailureModeSvc );
  
  DECLARE_ALGORITHM( CalXtalRecAlg );
  DECLARE_ALGORITHM( CalTupleAlg );

  DLL_DECL_TOOL( XtalDigiTool );
  DLL_DECL_TOOL( XtalRecTool );
  DLL_DECL_TOOL( CalTrigTool );
  DLL_DECL_TOOL( PrecalcCalibTool );
  
  DLL_DECL_TOOL( CalXtalRespRandom );
} 
