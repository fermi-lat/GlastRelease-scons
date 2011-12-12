// $Header$

/** @class XmlCalTholdCICnv

  Converter from xml to TCDS CalTholdCI class

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

#include "CalibData/Cal/CalTholdCI.h"

// Probably don't need the following:
// #include "CalibData/CalibTime.h"

#include "idents/CalXtalId.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

template <class TYPE> class CnvFactory;

class XmlCalTholdCICnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalTholdCICnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalTholdCICnv(ISvcLocator* svcs);

  virtual ~XmlCalTholdCICnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

private:
  /// Private utility which knows how to get the information out of a
  /// <tholdCI> element and make a CalibData::CalTholdCI with it
  CalibData::CalTholdCI* processRange(DOMElement* tholdCIElt);

};

// Begin implementation
//static CnvFactory<XmlCalTholdCICnv> s_factory;
//const  ICnvFactory& XmlCalTholdCICnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalTholdCICnv);



XmlCalTholdCICnv::XmlCalTholdCICnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_TholdCI) { 
}


const CLID& XmlCalTholdCICnv::objType() const {
  return CLID_Calib_CAL_TholdCI;
}

const CLID& XmlCalTholdCICnv::classID() {
  return CLID_Calib_CAL_TholdCI;
}

CalibData::CalTholdCI* 
XmlCalTholdCICnv::processRange(DOMElement* tholdCIElt) {
  using xmlBase::Dom;
  using CalibData::ValSig;
  using CalibData::CalTholdCI;
  using idents::CalXtalId;

  MsgStream log(msgSvc(), "XmlCalTholdCICnv" );
  
  // Could check here to make sure it really is a <tholdCI>
  
  ValSig* FLE = 0;
  ValSig* FHE = 0;
  ValSig* LAC = 0;
  // The "4" is for the four possible ranges
  std::vector<ValSig>* peds = new std::vector<ValSig>(4);
  std::vector<ValSig>* ULD  = new std::vector<ValSig>(4);
  
  try {
    FLE = processValSig(tholdCIElt, "FLEVal", "FLESig");
    FHE = processValSig(tholdCIElt, "FHEVal", "FHESig");
    LAC = processValSig(tholdCIElt, "LACVal", "LACSig");
  }
  catch (xmlBase::DomException ex1) {
    delete FLE; delete FHE; delete LAC;
    log << MSG::ERROR <<  ex1.getMsg() << endreq;
    return 0;
  }

  // Now for up to 4 tholdCIRange child elements
  DOMElement* rangeElt = Dom::getFirstChildElement(tholdCIElt);
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
      errMsg = "Bad range value for <tholdCI> child ";
      break;
    }
    if ((1 << ix) & mask) {
      problem = true;
      errMsg = "Duplicate range information in TholdCI calibration";
      break;
    }
    mask |= 1 << ix;

    // process
    try {
      ValSig* vs = processValSig(rangeElt, "pedVal", "pedSig");
      (*peds)[ix] = *vs;
      delete vs; vs = 0;
      vs = processValSig(rangeElt, "ULDVal", "ULDSig");
      (*ULD)[ix] = *vs;
      delete vs;
      rangeElt = Dom::getSiblingElement(rangeElt);
    }
    catch (xmlBase::DomException ex2) {
      problem = true;
      errMsg = ex2.getMsg();
      break;
    }
  }

  CalTholdCI *pNew = 0;
  if (problem) {
    log << MSG::ERROR << errMsg << endreq;
  }      
  else {
    pNew = new CalTholdCI(ULD, FLE, FHE, LAC, peds);
  }
  delete ULD, delete FLE; delete FHE; delete LAC; delete peds;

  return pNew;
}

// Create our specific object
StatusCode XmlCalTholdCICnv::i_createObj(const DOMElement* docElt, 
                                           DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::CalTholdCI;
  using CalibData::CalTholdCICol;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalTholdCICol* pObj = 
    new CalTholdCICol(nRow, nCol, nLayer, nXtal, nFace);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    CalTholdCI* pTholdCI = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pTholdCI);
    delete pTholdCI;
    rangeElt = findNextRange(rangeElt);
  }

  return StatusCode::SUCCESS;
}



