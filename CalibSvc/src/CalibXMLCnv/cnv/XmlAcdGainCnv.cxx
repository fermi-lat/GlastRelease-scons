// $Header$

#include <string>
#include "XmlAcdGainCnv.h"

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

#include "CalibData/Acd/AcdCalibGain.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

static CnvFactory<XmlAcdGainCnv> s_factory;
const  ICnvFactory& XmlAcdGainCnvFactory = s_factory;



XmlAcdGainCnv::XmlAcdGainCnv( ISvcLocator* svc) :
  XmlAcdBaseCnv(svc, CLID_Calib_ACD_ElecGain) { 
}


const CLID& XmlAcdGainCnv::objType() const {
  return CLID_Calib_ACD_ElecGain;
}

const CLID& XmlAcdGainCnv::classID() {
  return CLID_Calib_ACD_ElecGain;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <acdGain> element and make a CalibData::AcdGain with it
  CalibData::AcdGain* processRange(DOMElement* gainElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is an <acdGain>
    float gain, sig;
    try {
      gain = xmlBase::Dom::getDoubleAttribute(gainElt, "avg");
      sig = xmlBase::Dom::getDoubleAttribute(gainElt, "sig");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdGainCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AcdGain(gain, sig);
  }
}

// Create our specific object
StatusCode XmlAcdGainCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AcdGain;

  unsigned nFace, nRow, nCol, nPmt, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readAcdDimension(docElt, nFace, nRow, nCol, nPmt, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AcdCalibGain* pObj = 
    new CalibData::AcdCalibGain(nFace, nRow, nCol, nPmt, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    AcdGain* pGain = processRange(rangeElt);
    pObj->putRange(m_nFace, m_nRow, m_nCol, m_nPmt, m_nRange, 
                   pGain);
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
