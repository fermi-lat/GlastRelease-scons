// $Header$

/** @class XmlCalMevPerDacCnv

  Converter from xml to TCDS CalMevPerDac class

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

#include "CalibData/Cal/CalMevPerDac.h"

// Probably need the following:
#include "CalibData/Cal/Xpos.h"

#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

template <class TYPE> class CnvFactory;

class XmlCalMevPerDacCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalMevPerDacCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalMevPerDacCnv(ISvcLocator* svcs);

  virtual ~XmlCalMevPerDacCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

private:
  /// Private utility which knows how to get the information out of a
  /// <mevPerDac> element and make a CalibData::CalMevPerDac with it
  CalibData::CalMevPerDac* processRange(DOMElement* mevPerDacElt);

};

// Begin implementation
//static CnvFactory<XmlCalMevPerDacCnv> s_factory;
//const  ICnvFactory& XmlCalMevPerDacCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalMevPerDacCnv);



XmlCalMevPerDacCnv::XmlCalMevPerDacCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_MevPerDac) { 
}


const CLID& XmlCalMevPerDacCnv::objType() const {
  return CLID_Calib_CAL_MevPerDac;
}

const CLID& XmlCalMevPerDacCnv::classID() {
  return CLID_Calib_CAL_MevPerDac;
}

CalibData::CalMevPerDac* 
XmlCalMevPerDacCnv::processRange(DOMElement* mevPerDacElt) {
  using xmlBase::Dom;
  using CalibData::ValSig;
  using CalibData::CalMevPerDac;

  MsgStream log(msgSvc(), "XmlCalMevPerDacCnv" );
  
  // Could check here to make sure it really is a <mevPerDac>
  
  ValSig* big = 0;
  ValSig* small = 0;
  std::vector<ValSig>* bigSmallRatioN = 0;
  std::vector<ValSig>* bigSmallRatioP = 0;
  
  try {
    big = processValSig(mevPerDacElt, "bigVal", "bigSig");
    small = processValSig(mevPerDacElt, "smallVal", "smallSig");
    
    DOMElement* bigSmall = Dom::getFirstChildElement(mevPerDacElt);
    std::string face = Dom::getAttribute(bigSmall, "end");
    bool posEnd1 = (face == "POS");
    
    if (posEnd1) {
      bigSmallRatioP = processValSigs(bigSmall, 
                                      "bigSmallRatioVals",
                                      "bigSmallRatioSigs");
    }
    else {   // must be "NEG"
      bigSmallRatioN = processValSigs(bigSmall, 
                                      "bigSmallRatioVals",
                                      "bigSmallRatioSigs");
    }
    bigSmall = Dom::getSiblingElement(bigSmall);
    bool posEnd2 = (Dom::getAttribute(bigSmall, "end") == "POS");
    if (posEnd1 == posEnd2) {
      std::cerr << "From CalibSvc::XmlCalMevPerDacCnv::processRange  " 
                << "Illegal <mevPerDac> element " << std::endl;
      delete big; delete small; 
      if (bigSmallRatioN) delete bigSmallRatioN;
      if (bigSmallRatioP) delete bigSmallRatioP;
      return 0;
    }
    
    if (posEnd2) {
      bigSmallRatioP = processValSigs(bigSmall, 
                                      "bigSmallRatioVals",
                                      "bigSmallRatioSigs");
    }
    else {   // must be "NEG"
      bigSmallRatioN = processValSigs(bigSmall, 
                                      "bigSmallRatioVals",
                                      "bigSmallRatioSigs");
    }
  }
  catch (xmlBase::DomException ex) {
    log << MSG::ERROR  << ex.getMsg() << endreq;
    if (big) delete big;
    if (small) delete small;
    if (bigSmallRatioN) delete bigSmallRatioN; 
    if (bigSmallRatioP) delete bigSmallRatioP;
    throw ex;
  }
  
  CalMevPerDac *pNew = new CalMevPerDac(big, small, bigSmallRatioN,
                                        bigSmallRatioP);
  delete big; delete small;
  delete bigSmallRatioN; delete bigSmallRatioP;
  
  return pNew;
}

// Create our specific object
StatusCode XmlCalMevPerDacCnv::i_createObj(const DOMElement* docElt, 
                                           DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::CalMevPerDac;
  using CalibData::CalMevPerDacCol;
  using CalibData::Xpos;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange, nDacCol, nXpos;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange,
                  &nDacCol, &nXpos);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalMevPerDacCol* pObj = 
    new CalMevPerDacCol(nRow, nCol, nLayer, nXtal, nXpos);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    CalMevPerDac* pMevPerDac = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pMevPerDac);
    delete pMevPerDac;
    rangeElt = findNextRange(rangeElt);
  }

  // Also have to handle xpos if present
  DOMElement* xposElt = findXpos(docElt);

  if (xposElt != 0 ) {
    Xpos* pXpos = processXpos(xposElt);
    pObj->putXpos(pXpos);
    delete pXpos;
  }
  return StatusCode::SUCCESS;
}



