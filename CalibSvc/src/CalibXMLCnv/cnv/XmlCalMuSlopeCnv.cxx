// $Header$

#include <string>
#include "XmlCalMuSlopeCnv.h"

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

#include "CalibData/Cal/CalCalibMuSlope.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlCalMuSlopeCnv> s_factory;
//const  ICnvFactory& XmlCalMuSlopeCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalMuSlopeCnv);



XmlCalMuSlopeCnv::XmlCalMuSlopeCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_MuSlope) { 
}


const CLID& XmlCalMuSlopeCnv::objType() const {
  return CLID_Calib_CAL_MuSlope;
}

const CLID& XmlCalMuSlopeCnv::classID() {
  return CLID_Calib_CAL_MuSlope;
}

// Don't need to override in this case
/*
StatusCode XmlBaseCnv::i_processObj(DataObject*, // pObject,
                                    IOpaqueAddress*)  {  //pAddress
  return StatusCode::SUCCESS;
}
*/

namespace {
  /// Local utility which knows how to get the information out of a
  /// <calGain> element and make a CalibData::Gain with it
  CalibData::MuSlope* processRange(DOMElement* slopeElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is a <calMuSlope>
    float slope, error;

    try {
      slope = xmlBase::Dom::getDoubleAttribute(slopeElt, "slope");
      error = xmlBase::Dom::getDoubleAttribute(slopeElt, "error");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalMuySlopeCnv::processRange" 
                << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    /*
    std::string att = Dom::getAttribute(slopeElt, "slope");
    float slope = atof(att.c_str());
    att = Dom::getAttribute(slopeElt, "error");
    float error = atof(att.c_str());
    */
    return new CalibData::MuSlope(slope, error);
  }
}

// Create our specific object
StatusCode XmlCalMuSlopeCnv::i_createObj(const DOMElement* docElt, 
                                         DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::MuSlope;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibMuSlope* pObj = 
    new CalibData::CalCalibMuSlope(nRow, nCol, nLayer, nXtal, nFace, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    MuSlope* pMuSlope = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, m_nFace, 
                   pMuSlope);
    delete pMuSlope;
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
