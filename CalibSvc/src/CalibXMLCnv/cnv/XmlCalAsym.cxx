// $Header$

/** @class XmlCalAsymCnv

  Converter from xml to TCDS CalAsym class

  @author J. Bogart
*/
#include "XmlCalBaseCnv.h"
// #include <xercesc/dom/DOMElement.hpp>

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
#include "CalibData/Cal/Xpos.h"

#include "CalibData/Cal/CalAsym.h"

#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

template <class TYPE> class CnvFactory;

class XmlCalAsymCnv : public XmlCalBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<XmlCalAsymCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlCalAsymCnv(ISvcLocator* svcs);

  virtual ~XmlCalAsymCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

private:
  /// Utility which knows how to get the information out of an
  /// <asym> element and make a CalibData::CalAsym with it
  CalibData::CalAsym* processRange(DOMElement* mevPerDacElt);

};

// Begin implementation
//static CnvFactory<XmlCalAsymCnv> s_factory;
//const  ICnvFactory& XmlCalAsymCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalAsymCnv);



XmlCalAsymCnv::XmlCalAsymCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_Asym) { 
}


const CLID& XmlCalAsymCnv::objType() const {
  return CLID_Calib_CAL_Asym;
}

const CLID& XmlCalAsymCnv::classID() {
  return CLID_Calib_CAL_Asym;
}

CalibData::CalAsym* XmlCalAsymCnv::processRange(DOMElement* asymElt) {
  using xmlBase::Dom;
  using CalibData::ValSig;
  using CalibData::CalAsym;
  MsgStream log(msgSvc(), "XmlCalAsymCnv" );

  // Could check here to make sure it really is an <asym>

  std::vector<ValSig>* bigs = 0;
  std::vector<ValSig>* smalls = 0;
  std::vector<ValSig>* nSmallPBig = 0;
  std::vector<ValSig>* pSmallNBig = 0;

      
  try {
    bigs = processValSigs(asymElt, "bigVals", "bigSigs");
    smalls = processValSigs(asymElt, "smallVals", "smallSigs");
    nSmallPBig = processValSigs(asymElt, "NsmallPbigVals", "NsmallPbigSigs");
    pSmallNBig = processValSigs(asymElt, "PsmallNbigVals", "PsmallNbigSigs");

  }
  catch (xmlBase::DomException ex) {
    log << MSG::ERROR << "Error in calibration XML input" 
        << ex.getMsg() << endreq;
    if (bigs) delete bigs;
    if (smalls) delete smalls;
    if (nSmallPBig) delete nSmallPBig;
    if (pSmallNBig) delete pSmallNBig;
    throw ex;
  }

  CalAsym *pNew = new CalAsym(bigs, smalls, nSmallPBig, pSmallNBig);

  delete bigs; delete smalls;
  delete nSmallPBig;   delete pSmallNBig;
  
  return pNew;
}

// Create our specific object
StatusCode XmlCalAsymCnv::i_createObj(const DOMElement* docElt, 
                                           DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::CalAsym;
  using CalibData::CalAsymCol;
  using CalibData::Xpos;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange, nDacCol, nXpos;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange,
                  &nDacCol, &nXpos);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalAsymCol* pObj = 
    new CalAsymCol(nRow, nCol, nLayer, nXtal, nXpos);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    CalAsym* pAsym = processRange(rangeElt);
    pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                   m_nFace, pAsym);
    delete pAsym;
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



