// $Header$

#include <string>
#include "XmlCalGainCnv.h"

// #include <util/XMLUni.hpp>
// #include <util/XMLString.hpp>

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
#include "xml/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

static CnvFactory<XmlCalGainCnv> s_factory;
const  ICnvFactory& XmlCalGainCnvFactory = s_factory;



XmlCalGainCnv::XmlCalGainCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_CAL_ElecGain) { 
}


const CLID& XmlCalGainCnv::objType() const {
  return CLID_Calib_CAL_ElecGain;
}

const CLID& XmlCalGainCnv::classID() {
  return CLID_Calib_CAL_ElecGain;
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
  CalibData::Gain* processRange(DOM_Element gainElt) {
    using xml::Dom;

    // Could check here to make sure it really is a <calGain>

    std::string att = Dom::getAttribute(gainElt, "gain");
    float gain = atof(att.c_str());

    return new CalibData::Gain(gain);
  }
}

// Create our specific object
StatusCode XmlCalGainCnv::i_createObj(const DOM_Element& docElt, 
                                     DataObject*& refpObject)
{
  using xml::Dom;
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

  DOM_Element rangeElt = findFirstRange(docElt);

  while (rangeElt != DOM_Element() ) {
    Gain* pGain = processRange(rangeElt);
    pObj->putRange(nRow, nCol, nLayer, nXtal, nRange, nFace, pGain);
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
