//  $Header$
/**
   file XmlTkrSplitsCnv.cxx

   Conversion from persistent XML representation to calibration TDS
   for Tkr splits calibration
*/
#include <string>

#include "XmlTkrBaseCnv.h"

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

#include "CalibData/Tkr/TkrSplitsCalib.h"
#include "CalibData/CalibTime.h"
#include "CalibData/CalibModel.h"         // for definition of CLid

#include "xmlBase/Dom.h"


template <class TYPE> class CnvFactory;

class XmlTkrSplitsCnv : public XmlTkrBaseCnv {
  /// Friend needed for instantiation
  friend class CnvFactory<XmlTkrSplitsCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlTkrSplitsCnv(ISvcLocator* svcs);

  virtual ~XmlTkrSplitsCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);

};

//static CnvFactory<XmlTkrSplitsCnv> s_factory;
//const  ICnvFactory& XmlTkrSplitsCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlTkrSplitsCnv);

XmlTkrSplitsCnv::XmlTkrSplitsCnv( ISvcLocator* svc) :
  XmlTkrBaseCnv(svc, CLID_Calib_TKR_Splits) { 
}

const CLID& XmlTkrSplitsCnv::objType() const {
  return CLID_Calib_TKR_Splits;
}

const CLID& XmlTkrSplitsCnv::classID() {
  return CLID_Calib_TKR_Splits;
}


StatusCode XmlTkrSplitsCnv::i_createObj(const DOMElement* docElt, 
                                        DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::TkrSplit;

  unsigned nRow, nCol, nTray, nChip;
  StatusCode status = readDimension(docElt, nRow, nCol, nTray, nChip);
  if (!status.isSuccess()) return status;

  CalibData::TkrSplitsCalib* pObj =
    new CalibData::TkrSplitsCalib(nRow, nCol, nTray);

  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  std::vector<DOMElement*> splitsElts;

  Dom::getDescendantsByTagName(docElt, "split", splitsElts);
  
  unsigned nElt = splitsElts.size();

  for (unsigned iElt = 0; iElt < nElt; iElt++) {
    DOMElement* elt = splitsElts[iElt];
    try {
      unsigned iTowerRow = Dom::getIntAttribute(elt, "iBayRow");
      unsigned iTowerCol = Dom::getIntAttribute(elt, "iBayCol");
      unsigned iTray = Dom::getIntAttribute(elt, "iTray");
      unsigned low = Dom::getIntAttribute(elt, "low");
      unsigned high = Dom::getIntAttribute(elt, "high");
      std::string topString = Dom::getAttribute(elt, "top");
      bool top = (topString == "true");
      CalibData::TkrSplit split(low, high);
      pObj->putChannel(&split, iTowerRow, iTowerCol, iTray, top);
    }
    catch (xmlBase::DomException ex) {
      return StatusCode::FAILURE;
    }
  }
  return StatusCode::SUCCESS;
}
