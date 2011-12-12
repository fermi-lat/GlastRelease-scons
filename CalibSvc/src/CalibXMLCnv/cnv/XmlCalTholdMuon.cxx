// $Header$

/** @class XmlCalTholdMuonCnv

  Converter from xml to TCDS CalTholdMuon class

  @author J. Bogart
*/
#include "XmlCalBaseCnv.h"

// Following only needed for implementation, not class definition
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

#include "CalibData/Cal/CalTholdMuon.h"

// Probably don't need the following:
// #include "CalibData/CalibTime.h"

#include "idents/CalXtalId.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

template <class TYPE> class CnvFactory;

class XmlCalTholdMuonCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalTholdMuonCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalTholdMuonCnv(ISvcLocator* svcs);

  virtual ~XmlCalTholdMuonCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

private:
  /// Private utility which knows how to get the information out of a
  /// <tholdMuon> element and make a CalibData::CalTholdMuon with it
  CalibData::CalTholdMuon* processRange(DOMElement* tholdMuonElt);

};

// Begin implementation
//static CnvFactory<XmlCalTholdMuonCnv> s_factory;
//const  ICnvFactory& XmlCalTholdMuonCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalTholdMuonCnv);



XmlCalTholdMuonCnv::XmlCalTholdMuonCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_TholdMuon) { 
}


const CLID& XmlCalTholdMuonCnv::objType() const {
  return CLID_Calib_CAL_TholdMuon;
}

const CLID& XmlCalTholdMuonCnv::classID() {
  return CLID_Calib_CAL_TholdMuon;
}

CalibData::CalTholdMuon* 
XmlCalTholdMuonCnv::processRange(DOMElement* tholdMuonElt) {
  using xmlBase::Dom;
  using CalibData::ValSig;
  using CalibData::CalTholdMuon;
  using idents::CalXtalId;

  MsgStream log(msgSvc(), "XmlCalTholdMuonCnv" );
  
  // Could check here to make sure it really is a <tholdMuon>
  
  ValSig* FLE = 0;
  ValSig* FHE = 0;
  // The "4" is for the four possible ranges
  std::vector<ValSig>* peds = new std::vector<ValSig>(4);
  
  try {
    FLE = processValSig(tholdMuonElt, "FLEVal", "FLESig");
    FHE = processValSig(tholdMuonElt, "FHEVal", "FHESig");
  }
  catch (xmlBase::DomException ex1) {
    delete FLE; delete FHE;
    log << MSG::ERROR <<  ex1.getMsg() << endreq;
    return 0;
  }

  // Now for up to 4 tholdMuonRange child elements
  DOMElement* rangeElt = Dom::getFirstChildElement(tholdMuonElt);
  unsigned mask = 0;
  bool problem = false;
  std::string errMsg("");

  while (rangeElt != 0 ) {

    // check we don't already have this one
    std::string range = Dom::getAttribute(rangeElt, "range");
    unsigned ix;
    if (range == "LEX8") ix = CalXtalId::LEX8;
    else if (range == "LEX1" ) ix = CalXtalId::LEX1;
    else if (range == "HEX8" ) ix = CalXtalId::HEX8;
    else if (range == "HEX1" ) ix = CalXtalId::HEX1;
    else { 
      problem = true;
      errMsg = "Bad range value for <tholdMuon> child ";
      break;
    }
    if ((1 << ix) & mask) {
      problem = true;
      errMsg = "Duplicate range information in TholdMuon calibration";
      break;
    }
    mask |= 1 << ix;

    // process
    try {
      ValSig* ped = processValSig(rangeElt, "pedVal", "pedSig");
      (*peds)[ix] = *ped;
      delete ped;
      rangeElt = Dom::getSiblingElement(rangeElt);
    }
    catch (xmlBase::DomException ex2) {
      problem = true;
      errMsg = ex2.getMsg();
      break;
    }
  }

  CalTholdMuon *pNew = 0;
  if (problem) {
    log << MSG::ERROR << errMsg << endreq;
  }      
  else {
    pNew = new CalTholdMuon(FLE, FHE, peds);
  }
  delete FLE; delete FHE; delete peds;

  return pNew;
}

// Create our specific object
StatusCode XmlCalTholdMuonCnv::i_createObj(const DOMElement* docElt, 
                                           DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::CalTholdMuon;
  using CalibData::CalTholdMuonCol;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalTholdMuonCol* pObj = 
    new CalTholdMuonCol(nRow, nCol, nLayer, nXtal, nFace);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    CalTholdMuon* pTholdMuon = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pTholdMuon);
    delete pTholdMuon;
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}



