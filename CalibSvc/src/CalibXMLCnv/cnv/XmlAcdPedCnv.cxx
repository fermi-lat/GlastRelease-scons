// $Header$

#include <string>
#include "XmlAcdPedCnv.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibSvc/ICalibXmlSvc.h"
#include "CalibSvc/ICalibMetaCnvSvc.h"

#include "CalibData/Acd/AcdCalibPed.h"
#include "CalibData/CalibTime.h"
#include "xml/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

static CnvFactory<XmlAcdPedCnv> s_factory;
const  ICnvFactory& XmlAcdPedCnvFactory = s_factory;



XmlAcdPedCnv::XmlAcdPedCnv( ISvcLocator* svc) :
  XmlAcdBaseCnv(svc, CLID_Calib_ACD_Ped) { 
}


const CLID& XmlAcdPedCnv::objType() const {
  return CLID_Calib_ACD_Ped;
}

const CLID& XmlAcdPedCnv::classID() {
  return CLID_Calib_ACD_Ped;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <acdPed> element and make a CalibData::AcdPed with it
  CalibData::AcdPed* processRange(DOM_Element pedElt) {
    using xml::Dom;

    // Could check here to make sure it really is an <acdPed>
    float ped, sig;
    try {
      ped = xml::Dom::getDoubleAttribute(pedElt, "avg");
      sig = xml::Dom::getDoubleAttribute(pedElt, "sig");
    }
    catch (xml::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdPedCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AcdPed(ped, sig);
  }
}

// Create our specific object
StatusCode XmlAcdPedCnv::i_createObj(const DOM_Element& docElt, 
                                     DataObject*& refpObject)
{
  using xml::Dom;
  using CalibData::AcdPed;

  unsigned nFace, nRow, nCol, nPmt, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readAcdDimension(docElt, nFace, nRow, nCol, nPmt, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AcdCalibPed* pObj = 
    new CalibData::AcdCalibPed(nFace, nRow, nCol, nPmt, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOM_Element rangeElt = findFirstRange(docElt);

  while (rangeElt != DOM_Element() ) {
    AcdPed* pPed = processRange(rangeElt);
    pObj->putRange(m_nFace, m_nRow, m_nCol, m_nPmt, m_nRange, 
                   pPed);
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
