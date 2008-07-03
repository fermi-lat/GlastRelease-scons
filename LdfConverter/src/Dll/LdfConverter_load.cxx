/** 
* @file LdfConverter_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(LdfConverter) {
  DECLARE_SERVICE(LdfCnvSvc);
  DECLARE_SERVICE( LdfEventSelector );
  DECLARE_CONVERTER(LdfEventSummaryCnv);
  DECLARE_CONVERTER(LdfEventCnv);
  DECLARE_CONVERTER(LdfDigiCnv);
  DECLARE_CONVERTER(LdfTkrDigiCnv);
  DECLARE_CONVERTER(LdfCalDigiCnv);
  DECLARE_CONVERTER(LdfAcdDigiCnv);
  DECLARE_CONVERTER(LdfDiagnosticCnv);
  DECLARE_CONVERTER(LdfTimeCnv);
  DECLARE_CONVERTER(LdfGemCnv);
  DECLARE_CONVERTER( McEventCnv );
  DECLARE_CONVERTER( MetaEventCnv );
  DECLARE_CONVERTER( CcsdsCnv );
  DECLARE_CONVERTER( AncillaryEventCnv );
  DECLARE_CONVERTER( AdfEventCnv );

  DECLARE_CONVERTER( FilterCnv ); 
  DECLARE_CONVERTER( EbfCnv ); 

}
  
