/** 
* @file CalibSvc_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(CalibSvc) {
  // Useful for now to fake event time  
  DECLARE_ALGORITHM(CalibEvtClock);

  DECLARE_SERVICE(CalibDataSvc);


  DECLARE_SERVICE(CalibMySQLCnvSvc);

  DECLARE_SERVICE(CalibXmlCnvSvc);

  DECLARE_CONVERTER(XmlTest1Cnv);
  DECLARE_CONVERTER(XmlBadStripsCnv);
  DECLARE_CONVERTER(XmlCalPedCnv);
  DECLARE_CONVERTER(XmlCalGainCnv);
  DECLARE_CONVERTER(XmlCalMuSlopeCnv);
  DECLARE_CONVERTER(XmlCalLightAttCnv);
  DECLARE_CONVERTER(XmlCalLightAsymCnv);
  DECLARE_CONVERTER(XmlCalIntNonlinCnv);

  DECLARE_CONVERTER(XmlCalMevPerDacCnv);
  DECLARE_CONVERTER(XmlCalAsymCnv);
  DECLARE_CONVERTER(XmlCalTholdMuonCnv);
  DECLARE_CONVERTER(XmlCalTholdCICnv);

  DECLARE_CONVERTER(XmlAcdPedCnv);
  DECLARE_CONVERTER(XmlAcdGainCnv);

  DECLARE_CONVERTER(XmlTkrSplitsCnv);
  DECLARE_CONVERTER(XmlTkrTowerAlignCnv);
  DECLARE_CONVERTER(XmlTkrInternalAlignCnv);

  DECLARE_CONVERTER(XmlSAABoundaryCnv);
  DECLARE_CONVERTER(XmlLATAlignmentCnv);

  DECLARE_SERVICE(CalibRootCnvSvc);
  DECLARE_CONVERTER(RootCalGainCnv);
  DECLARE_CONVERTER(RootTkrTotCnv);
  DECLARE_CONVERTER(RootTkrScaleCnv);
} 
