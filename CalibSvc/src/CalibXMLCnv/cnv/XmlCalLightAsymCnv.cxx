// $Header$

#include <string>
#include "XmlCalLightAsymCnv.h"

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

#include "CalibData/Cal/CalCalibLightAsym.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlCalLightAsymCnv> s_factory;
//const  ICnvFactory& XmlCalLightAsymCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalLightAsymCnv);



XmlCalLightAsymCnv::XmlCalLightAsymCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_LightAsym) { 
}


const CLID& XmlCalLightAsymCnv::objType() const {
  return CLID_Calib_CAL_LightAsym;
}

const CLID& XmlCalLightAsymCnv::classID() {
  return CLID_Calib_CAL_LightAsym;
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
  /// <calLightAsym> element and make a CalibData::LightAsym with it
  CalibData::LightAsym* processRange(DOMElement* lightAsymElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is a <lightAsym>

    float err;
    std::vector<float> vals;
    try {
      xmlBase::Dom::getFloatsAttribute(lightAsymElt, "values", vals);
      err = xmlBase::Dom::getDoubleAttribute(lightAsymElt, "error");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalLightAsymCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::LightAsym(&vals, err);
  }
}

// Create our specific object
StatusCode XmlCalLightAsymCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::LightAsym;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibLightAsym* pObj = 
    new CalibData::CalCalibLightAsym(nRow, nCol, nLayer, nXtal, nFace, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    LightAsym* pLightAsym = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pLightAsym);
    delete pLightAsym;
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
