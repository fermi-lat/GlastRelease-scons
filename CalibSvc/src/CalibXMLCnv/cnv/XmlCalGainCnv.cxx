// $Header$

#include <string>
#include "XmlCalGainCnv.h"

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

#include "CalibData/Cal/CalCalibGain.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlCalGainCnv> s_factory;
//const  ICnvFactory& XmlCalGainCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalGainCnv);



XmlCalGainCnv::XmlCalGainCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_ElecGain) { 
}


const CLID& XmlCalGainCnv::objType() const {
  return CLID_Calib_CAL_ElecGain;
}

const CLID& XmlCalGainCnv::classID() {
  return CLID_Calib_CAL_ElecGain;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <calGain> element and make a CalibData::Gain with it
  CalibData::Gain* processRange(DOMElement* gainElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is a <calGain>
    float gain, sig;
    try {
      gain = xmlBase::Dom::getDoubleAttribute(gainElt, "avg");
      sig = xmlBase::Dom::getDoubleAttribute(gainElt, "sig");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalGainCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::Gain(gain, sig);
  }
}

// Create our specific object
StatusCode XmlCalGainCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::Gain;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibGain* pObj = 
    new CalibData::CalCalibGain(nRow, nCol, nLayer, nXtal, nFace, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    Gain* pGain = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                             m_nFace, pGain);
    delete pGain;
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
