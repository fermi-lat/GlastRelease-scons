/** 
* @file CalibSvc_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(CalibSvc) {
  DECLARE_SERVICE(CalibDataSvc);


  DECLARE_SERVICE(CalibMySQLCnvSvc);
  DECLARE_SERVICE(CalibXmlCnvSvc);

  DECLARE_CONVERTER(XmlBaseCnv);  

  // Following don't exist yet.
  //  DECLARE_SERVICE(CalibRootCnvSvc);
  // DECLARE_CONVERTER(XmlTest1Cnv);


} 
