// $Header$

#include <string>

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

#include "CalibData/Anc/AncCalibQdcPed.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAncQdcPedCnv

  Converter from xml to TCDS ANC pedestals class

  @author J. Bogart
*/
#include "XmlAncBaseCnv.h"

//class XmlAncQdcPedCnv;

// template <class TYPE> class CnvFactory;
//static CnvFactory<XmlAncQdcPedCnv> s_factory;
//const  ICnvFactory& XmlAncQdcPedCnvFactory = s_factory;
//DECLARE_CONVERTER_FACTORY(XmlAncQdcPedCnv);

class XmlAncQdcPedCnv : public XmlAncBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAncQdcPedCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAncQdcPedCnv(ISvcLocator* svcs);

  virtual ~XmlAncQdcPedCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

DECLARE_CONVERTER_FACTORY(XmlAncQdcPedCnv);

//

XmlAncQdcPedCnv::XmlAncQdcPedCnv( ISvcLocator* svc) :
  XmlAncBaseCnv(svc, CLID_Calib_ANC_QdcPed) { 
}


const CLID& XmlAncQdcPedCnv::objType() const {
  return CLID_Calib_ANC_QdcPed;
}

const CLID& XmlAncQdcPedCnv::classID() {
  return CLID_Calib_ANC_QdcPed;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <qdcPed> element and make a CalibData::AncQdcPed with it
  CalibData::AncQdcPed* processChan(DOMElement* chanElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <chan>
    DOMElement* pedElt = xmlBase::Dom::getFirstChildElement(chanElt);

    // Could check here to make sure it really is an <ancPed>
    float val, rms;
    unsigned   isBad;
    std::string dev;
    try {
      dev = xmlBase::Dom::getAttribute(chanElt, "device");
      val = xmlBase::Dom::getDoubleAttribute(pedElt, "ped");
      rms = xmlBase::Dom::getDoubleAttribute(pedElt, "pedRms");
      isBad = xmlBase::Dom::getIntAttribute(pedElt, "isBad");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAncQdcPedCnv::processChan" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AncQdcPed(val, rms, isBad, dev);
  }
}

// Create our specific object
StatusCode XmlAncQdcPedCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AncQdcPed;

  unsigned nMod, nLay, nChan;
  // need dimensions to call the constructor
  StatusCode status = 
    readAncDimension(docElt, nMod, nLay, nChan);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AncCalibQdcPed* pObj = 
    new CalibData::AncCalibQdcPed(nMod, nChan);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* chanElt = findFirstChan(docElt);

  while (chanElt != 0 ) {
    AncQdcPed* pPed = processChan(chanElt);
    pObj->putChan(m_iMod, m_iLay, m_iChan, pPed);
    chanElt = findNextChan(chanElt);
  }

  return StatusCode::SUCCESS;
}
