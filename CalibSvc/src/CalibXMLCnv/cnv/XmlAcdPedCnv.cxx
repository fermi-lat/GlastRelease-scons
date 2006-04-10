// $Header$

#include <string>
//#include "XmlAcdPedCnv.h"

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
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAcdPedCnv

  Converter from xml to TCDS ACD pedestals class

  @author J. Bogart
*/
#include "XmlAcdBaseCnv.h"

class XmlAcdPedCnv;

// template <class TYPE> class CnvFactory;
static CnvFactory<XmlAcdPedCnv> s_factory;
const  ICnvFactory& XmlAcdPedCnvFactory = s_factory;

class XmlAcdPedCnv : public XmlAcdBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAcdPedCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAcdPedCnv(ISvcLocator* svcs);

  virtual ~XmlAcdPedCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};


//

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
  CalibData::AcdPed* processPmt(DOMElement* pmtElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <pmt>
    DOMElement* pedElt = xmlBase::Dom::getFirstChildElement(pmtElt);

    // Could check here to make sure it really is an <acdPed>
    float mean, width;
    unsigned   status;
    try {
      mean = xmlBase::Dom::getDoubleAttribute(pedElt, "mean");
      width = xmlBase::Dom::getDoubleAttribute(pedElt, "width");
      status = xmlBase::Dom::getIntAttribute(pedElt, "status");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdPedCnv::processPmt" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AcdPed(mean, width, status);
  }
}

// Create our specific object
StatusCode XmlAcdPedCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AcdPed;

  unsigned nFace, nRow, nCol, nPmt, nDet, nNA;
  // need dimensions to call the constructor
  StatusCode status = 
    readAcdDimension(docElt, nDet, nFace, nRow, nCol, nNA, nPmt);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AcdCalibPed* pObj = 
    new CalibData::AcdCalibPed(nFace, nRow, nCol, nNA, nPmt);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* pmtElt = findFirstPmt(docElt);

  while (pmtElt != 0 ) {
    AcdPed* pPed = processPmt(pmtElt);
    pObj->putPmt(m_id, m_nPmt, pPed);
    pmtElt = findNextPmt(pmtElt);
  }

  return StatusCode::SUCCESS;
}
