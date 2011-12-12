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

#include "CalibData/Anc/AncCalibTaggerPed.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAncTaggerPedCnv

  Converter from xml to TCDS ANC pedestals class

  @author J. Bogart
*/
#include "XmlAncBaseCnv.h"

//class XmlAncTaggerPedCnv;

// template <class TYPE> class CnvFactory;
//static CnvFactory<XmlAncTaggerPedCnv> s_factory;
//const  ICnvFactory& XmlAncTaggerPedCnvFactory = s_factory;

class XmlAncTaggerPedCnv : public XmlAncBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAncTaggerPedCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAncTaggerPedCnv(ISvcLocator* svcs);

  virtual ~XmlAncTaggerPedCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

DECLARE_CONVERTER_FACTORY(XmlAncTaggerPedCnv);

//

XmlAncTaggerPedCnv::XmlAncTaggerPedCnv( ISvcLocator* svc) :
  XmlAncBaseCnv(svc, CLID_Calib_ANC_TaggerPed) { 
}


const CLID& XmlAncTaggerPedCnv::objType() const {
  return CLID_Calib_ANC_TaggerPed;
}

const CLID& XmlAncTaggerPedCnv::classID() {
  return CLID_Calib_ANC_TaggerPed;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <taggerPed> element and make a CalibData::AncTaggerPed with it
  CalibData::AncTaggerPed* processChan(DOMElement* chanElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <chan>
    DOMElement* pedElt = xmlBase::Dom::getFirstChildElement(chanElt);

    // Could check here to make sure it really is an <ancPed>
    float val, rNoise, sNoise;
    unsigned   isBad;
    try {
      val = xmlBase::Dom::getDoubleAttribute(pedElt, "ped");
      rNoise = xmlBase::Dom::getDoubleAttribute(pedElt, "rNoise");
      sNoise = xmlBase::Dom::getDoubleAttribute(pedElt, "sNoise");
      isBad = xmlBase::Dom::getIntAttribute(pedElt, "isBad");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAncTaggerPedCnv::processChan" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AncTaggerPed(val, rNoise, sNoise, isBad);
  }
}

// Create our specific object
StatusCode XmlAncTaggerPedCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AncTaggerPed;

  unsigned nMod, nLay, nChan;
  // need dimensions to call the constructor
  StatusCode status = 
    readAncDimension(docElt, nMod, nLay, nChan);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AncCalibTaggerPed* pObj = 
    new CalibData::AncCalibTaggerPed(nMod, nLay, nChan);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* chanElt = findFirstChan(docElt);

  while (chanElt != 0 ) {
    AncTaggerPed* pPed = processChan(chanElt);
    pObj->putChan(m_iMod, m_iLay, m_iChan, pPed);
    chanElt = findNextChan(chanElt);
  }

  return StatusCode::SUCCESS;
}
