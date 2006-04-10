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

#include "CalibData/Acd/AcdCalibCno.h"
#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"



/** @class XmlAcdCnoCnv

  Converter from xml to TCDS ACD cno threshold class

  @author J. Bogart
*/
#include "XmlAcdBaseCnv.h"

class XmlAcdCnoCnv;

// template <class TYPE> class CnvFactory;
static CnvFactory<XmlAcdCnoCnv> s_factory;
const  ICnvFactory& XmlAcdCnoCnvFactory = s_factory;

class XmlAcdCnoCnv : public XmlAcdBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlAcdCnoCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlAcdCnoCnv(ISvcLocator* svcs);

  virtual ~XmlAcdCnoCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

XmlAcdCnoCnv::XmlAcdCnoCnv( ISvcLocator* svc) :
  XmlAcdBaseCnv(svc, CLID_Calib_ACD_ThreshHigh) { 
}


const CLID& XmlAcdCnoCnv::objType() const {
  return CLID_Calib_ACD_ThreshHigh;
}

const CLID& XmlAcdCnoCnv::classID() {
  return CLID_Calib_ACD_ThreshHigh;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <acdCno> element and make a CalibData::AcdCno with it
  CalibData::AcdCno* processPmt(DOMElement* pmtElt) {
    using xmlBase::Dom;

    // Element we're interested in is child of <pmt>
    DOMElement* cnoElt = xmlBase::Dom::getFirstChildElement(pmtElt);

    // Could check here to make sure it really is an <acdCno>
    float cno, width;
    unsigned   status;
    try {
      cno = xmlBase::Dom::getDoubleAttribute(cnoElt, "cno");
      width = xmlBase::Dom::getDoubleAttribute(cnoElt, "width");
      status = xmlBase::Dom::getIntAttribute(cnoElt, "status");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdCnoCnv::processPmt" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    return new CalibData::AcdCno(cno, width, status);
  }
}

// Create our specific object
StatusCode XmlAcdCnoCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::AcdCno;

  unsigned nFace, nRow, nCol, nPmt, nDet, nNA;
  // need dimensions to call the constructor
  StatusCode status = 
    readAcdDimension(docElt, nDet, nFace, nRow, nCol, nNA, nPmt);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::AcdCalibCno* pObj = 
    new CalibData::AcdCalibCno(nFace, nRow, nCol, nNA, nPmt);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* pmtElt = findFirstPmt(docElt);

  while (pmtElt != 0 ) {
    AcdCno* pCno = processPmt(pmtElt);
    pObj->putPmt(m_id, m_nPmt, pCno);
    pmtElt = findNextPmt(pmtElt);
  }

  return StatusCode::SUCCESS;
}
