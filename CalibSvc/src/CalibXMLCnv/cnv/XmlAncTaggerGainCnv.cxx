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

#include "CalibData/Anc/AncCalibTaggerGain.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAncTaggerGainCnv

  Converter from xml to TCDS ANC pedestals class

  @author J. Bogart
*/
#include "XmlAncBaseCnv.h"

//class XmlAncTaggerGainCnv;

// template <class TYPE> class CnvFactory;
//static CnvFactory<XmlAncTaggerGainCnv> s_factory;
//const  ICnvFactory& XmlAncTaggerGainCnvFactory = s_factory;

class XmlAncTaggerGainCnv : public XmlAncBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAncTaggerGainCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAncTaggerGainCnv(ISvcLocator* svcs);

  virtual ~XmlAncTaggerGainCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

DECLARE_CONVERTER_FACTORY(XmlAncTaggerGainCnv);

//

XmlAncTaggerGainCnv::XmlAncTaggerGainCnv( ISvcLocator* svc) :
  XmlAncBaseCnv(svc, CLID_Calib_ANC_TaggerGain) { 
}

const CLID& XmlAncTaggerGainCnv::objType() const {
  return CLID_Calib_ANC_TaggerGain;
}

const CLID& XmlAncTaggerGainCnv::classID() {
  return CLID_Calib_ANC_TaggerGain;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <taggerGain> element and make a CalibData::AncTaggerGain with it
  CalibData::AncTaggerGain* processChan(DOMElement* chanElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <chan>
    DOMElement* gainElt = xmlBase::Dom::getFirstChildElement(chanElt);

    // Could check here to make sure it really is an <ancGain>
    float val;
    unsigned   isBad;
    try {
      val = xmlBase::Dom::getDoubleAttribute(gainElt, "gain");
      isBad = xmlBase::Dom::getIntAttribute(gainElt, "isBad");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAncTaggerGainCnv::processChan" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AncTaggerGain(val, isBad);
  }
}

// Create our specific object
StatusCode XmlAncTaggerGainCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AncTaggerGain;

  unsigned nMod, nLay, nChan;
  // need dimensions to call the constructor
  StatusCode status = 
    readAncDimension(docElt, nMod, nLay, nChan);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AncCalibTaggerGain* pObj = 
    new CalibData::AncCalibTaggerGain(nMod, nLay, nChan);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* chanElt = findFirstChan(docElt);

  while (chanElt != 0 ) {
    AncTaggerGain* pGain = processChan(chanElt);
    pObj->putChan(m_iMod, m_iLay, m_iChan, pGain);
    chanElt = findNextChan(chanElt);
  }

  return StatusCode::SUCCESS;
}
