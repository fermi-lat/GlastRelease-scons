//====================================================================
//  DbCnv_load.cpp
//--------------------------------------------------------------------
//
//  Package    : DbCnv
//
//  Description: Implementation of <Package>_load routine. This routine 
//               is needed for forcing the linker to load all the components
//               of the library.. 
//
//====================================================================

#include "Gaudi/Interfaces/ICnvFactory.h"
#include "Gaudi/Interfaces/ISvcFactory.h"
#include "Gaudi/Interfaces/IAlgFactory.h"

#define DLL_DECL_OBJECT(x)     extern const IFactory&    x##Factory; x##Factory.addRef();
#define DLL_DECL_SERVICE(x)    extern const ISvcFactory& x##Factory; x##Factory.addRef();
#define DLL_DECL_CONVERTER(x)  extern const ICnvFactory& DbOOMs##x##Factory; DbOOMs##x##Factory.addRef();
#define DLL_DECL_ALGORITHM(x)  extern const IAlgFactory& x##Factory; x##Factory.addRef();

void DbConverters_load()     {

  // Declare all the services
  DLL_DECL_SERVICE( DbEventCnvSvc );
  DLL_DECL_SERVICE( DbEvtSelector );
  
  // Now declare all availible converters
  DLL_DECL_CONVERTER( EventCnv );
  DLL_DECL_CONVERTER( MCEventCnv );
  DLL_DECL_CONVERTER( MCACDHitVectorCnv );
   //DLL_DECL_CONVERTER( MCTKRHitVectorCnv );
  DLL_DECL_CONVERTER( MCSiLayerVectorCnv );
  DLL_DECL_CONVERTER( MCCalorimeterHitVectorCnv );
  //DLL_DECL_CONVERTER( AxPartCandidateVectorCnv );
  //DLL_DECL_CONVERTER( ExportConverter );
  //DLL_DECL_CONVERTER( ExportVectorCnv );
  //DLL_DECL_CONVERTER( ExportListCnv );
} 

extern "C" void DbConverters_loadRef() {
  DbConverters_load();
}
