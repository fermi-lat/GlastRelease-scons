// $Header$

#include <string>
#include "XmlCalPedCnv.h"

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

#include "CalibData/Cal/CalCalibPed.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlCalPedCnv> s_factory;
//const  ICnvFactory& XmlCalPedCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalPedCnv);



XmlCalPedCnv::XmlCalPedCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_Ped) { 
}


const CLID& XmlCalPedCnv::objType() const {
  return CLID_Calib_CAL_Ped;
}

const CLID& XmlCalPedCnv::classID() {
  return CLID_Calib_CAL_Ped;
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
  /// <calPed> element and make a CalibData::Ped with it
  CalibData::Ped* processRange(DOMElement* pedElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is a <calPed>
    float avg, sig, cos;
    try {
      avg = xmlBase::Dom::getDoubleAttribute(pedElt, "avg");
      sig = xmlBase::Dom::getDoubleAttribute(pedElt, "sig");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalPedCnv::processRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    // backward compatibility
    try {
      cos = xmlBase::Dom::getDoubleAttribute(pedElt, "cos");
    }
    catch (xmlBase::NullNode nullEx) {
      cos=2.0;
    }
    catch (xmlBase::DomException ex) {
        std::cerr << "From CalibSvc::XmlCalPedCnv::processRange" << std::endl;
        std::cerr << ex.getMsg() << std::endl;
        throw ex;
        //      }
    }
    return new CalibData::Ped(avg, sig, cos);
  }
}

// Create our specific object
StatusCode XmlCalPedCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::Ped;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibPed* pObj = 
    new CalibData::CalCalibPed(nRow, nCol, nLayer, nXtal, nFace, nRange);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    Ped* ped = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, 
                     m_nRange, m_nFace, ped);
    delete ped;

    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}
