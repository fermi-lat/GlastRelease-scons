// $Header$

#include <string>
#include "XmlCalIntNonlinCnv.h"

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

#include "CalibData/Cal/CalCalibIntNonlin.h"

// Probably need the following:
#include "CalibData/DacCol.h"

#include "CalibData/CalibTime.h"
#include "xml/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

static CnvFactory<XmlCalIntNonlinCnv> s_factory;
const  ICnvFactory& XmlCalIntNonlinCnvFactory = s_factory;



XmlCalIntNonlinCnv::XmlCalIntNonlinCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_IntNonlin) { 
}


const CLID& XmlCalIntNonlinCnv::objType() const {
  return CLID_Calib_CAL_IntNonlin;
}

const CLID& XmlCalIntNonlinCnv::classID() {
  return CLID_Calib_CAL_IntNonlin;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <calLightAsym> element and make a CalibData::IntNonlin with it
  CalibData::IntNonlin* processRange(DOM_Element intNonlinElt) {
    using xml::Dom;

    // Could check here to make sure it really is a <intNonlin>

    float err;
    std::vector<float> vals;
    try {
      xml::Dom::getFloatsAttribute(intNonlinElt, "values", vals);
      err = xml::Dom::getDoubleAttribute(intNonlinElt, "error");
    }
    catch (xml::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalIntNonlinCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::IntNonlin(&vals, err);
  }
}

// Create our specific object
StatusCode XmlCalIntNonlinCnv::i_createObj(const DOM_Element& docElt, 
                                     DataObject*& refpObject)
{
  using xml::Dom;
  using CalibData::IntNonlin;
  using CalibData::DacCol;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange, nDacCol;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange,
                  &nDacCol);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibIntNonlin* pObj = 
    new CalibData::CalCalibIntNonlin(nRow, nCol, nLayer, nXtal, nFace, nRange,
                                     nDacCol);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOM_Element rangeElt = findFirstRange(docElt);

  while (rangeElt != DOM_Element() ) {
    IntNonlin* pIntNonlin = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pIntNonlin);
    rangeElt = findNextRange(rangeElt);
  }

  // Also have to handle dac collections        <<<<-------
  DOM_Element dacColElt = findFirstDacCol(docElt);

  while (dacColElt != DOM_Element() ) {
    unsigned range;
    DacCol* pDacCol = processDacCol(dacColElt, &range);
    pObj->putDacCol(range, pDacCol);
    dacColElt = findNextDacCol(dacColElt);
  }
  return StatusCode::SUCCESS;
}


