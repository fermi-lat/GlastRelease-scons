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

#include "CalibData/Acd/AcdCalibVeto.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAcdVetoCnv

  Converter from xml to TCDS ACD veto threshold class

  @author J. Bogart
*/
#include "XmlAcdBaseCnv.h"

class XmlAcdVetoCnv;

// template <class TYPE> class CnvFactory;
static CnvFactory<XmlAcdVetoCnv> s_factory;
const  ICnvFactory& XmlAcdVetoCnvFactory = s_factory;

class XmlAcdVetoCnv : public XmlAcdBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAcdVetoCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAcdVetoCnv(ISvcLocator* svcs);

  virtual ~XmlAcdVetoCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

XmlAcdVetoCnv::XmlAcdVetoCnv( ISvcLocator* svc) :
  XmlAcdBaseCnv(svc, CLID_Calib_ACD_ThreshVeto) { 
}


const CLID& XmlAcdVetoCnv::objType() const {
  return CLID_Calib_ACD_ThreshVeto;
}

const CLID& XmlAcdVetoCnv::classID() {
  return CLID_Calib_ACD_ThreshVeto;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <acdVeto> element and make a CalibData::AcdVeto with it
  CalibData::AcdVeto* processPmt(DOMElement* pmtElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <pmt>
    DOMElement* vetoElt = xmlBase::Dom::getFirstChildElement(pmtElt);

    // Could check here to make sure it really is an <acdVeto>
    float veto, width;
    unsigned   status;
    try {
      veto = xmlBase::Dom::getDoubleAttribute(vetoElt, "veto");
      width = xmlBase::Dom::getDoubleAttribute(vetoElt, "width");
      status = xmlBase::Dom::getIntAttribute(vetoElt, "status");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdVetoCnv::processPmt" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AcdVeto(veto, width, status);
  }
}

// Create our specific object
StatusCode XmlAcdVetoCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AcdVeto;

  unsigned nFace, nRow, nCol, nPmt, nDet, nNA;
  // need dimensions to call the constructor
  StatusCode status = 
    readAcdDimension(docElt, nDet, nFace, nRow, nCol, nNA, nPmt);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AcdCalibVeto* pObj = 
    new CalibData::AcdCalibVeto(nFace, nRow, nCol, nNA, nPmt);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* pmtElt = findFirstPmt(docElt);

  while (pmtElt != 0 ) {
    AcdVeto* pVeto = processPmt(pmtElt);
    pObj->putPmt(m_id, m_nPmt, pVeto);
    pmtElt = findNextPmt(pmtElt);
  }

  return StatusCode::SUCCESS;
}
