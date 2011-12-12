//  $Header$
/**
   file XmlTkrAlignCnv.cxx

   Conversion from persistent XML representation to calibration TDS
   for Tkr alignment
*/
#include "XmlBaseCnv.h"

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

#include "CalibData/Tkr/TkrTowerAlignCalib.h"
#include "CalibData/Tkr/TkrInternalAlignCalib.h"

#include "CalibData/CalibTime.h"
#include "CalibData/CalibModel.h"         // for definition of CLid

#include "xmlBase/Dom.h"
#include "facilities/Util.h"

#define MAX_TOWER_ID 15

namespace {
  using CLHEP::Hep3Vector;
  using facilities::Util;
  using xmlBase::Dom;

  /**
     Input string @a att should be of form  "f1 f2 f3" where f1, f2 and f3 are
     valid floating point numbers.  Extract and store values in new
     Hep3Vector @a v.  Caller is responsible for deleting @a v
   */
  StatusCode stringTo3Vector(const std::string& att, Hep3Vector*& v) {
    std::vector<std::string> toks;
    Util::stringTokenize(att, " ", toks);

    if (toks.size() != 3)  return StatusCode::FAILURE;
    try {
      v = new Hep3Vector(Util::stringToDouble(toks[0]),
                         Util::stringToDouble(toks[1]),
                         Util::stringToDouble(toks[2]));
    }
    catch (std::exception ex) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
  /**
     Given DOM element @v elt  parse
     and store attributes id, disp and rot. Return failure if id < 0
     or if disp or rot are not properly formed.
   */
  StatusCode parseAlignElt(DOMElement* elt,  unsigned & id,
                           Hep3Vector*& pDisp, Hep3Vector*& pRot) {
    int intAtt;
    std::string dispAtt, rotAtt;
    try {
      dispAtt = Dom::getAttribute(elt, "disp");
      rotAtt = Dom::getAttribute(elt, "rot");
      intAtt = Dom::getIntAttribute(elt, "id");
    }
    catch (std::exception ex) {
      return StatusCode::FAILURE;
    }
    if (intAtt < 0)  return StatusCode::FAILURE;
    id = (unsigned) intAtt;
    pDisp = pRot = 0;
    StatusCode ret = stringTo3Vector(dispAtt, pDisp);
    if (ret == StatusCode::SUCCESS) {
      ret = stringTo3Vector(rotAtt, pRot);
    }
    return ret;
  }
}

template <class TYPE> class CnvFactory;

class XmlTkrTowerAlignCnv : public XmlBaseCnv {
  /// Friend needed for instantiation
  friend class CnvFactory<XmlTkrTowerAlignCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlTkrTowerAlignCnv(ISvcLocator* svcs);

  virtual ~XmlTkrTowerAlignCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);
};

//static CnvFactory<XmlTkrTowerAlignCnv> s_towerFactory;
//const  ICnvFactory& XmlTkrTowerAlignCnvFactory = s_towerFactory;
DECLARE_CONVERTER_FACTORY(XmlTkrTowerAlignCnv);

XmlTkrTowerAlignCnv::XmlTkrTowerAlignCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_TKR_TowerAlign) { }

const CLID& XmlTkrTowerAlignCnv::objType() const {
  return CLID_Calib_TKR_TowerAlign;
}

const CLID& XmlTkrTowerAlignCnv::classID() {
  return CLID_Calib_TKR_TowerAlign;
}

StatusCode XmlTkrTowerAlignCnv::i_createObj(const DOMElement* docElt, 
                                            DataObject*& refpObject)  {
  using xmlBase::Dom;
  using CalibData::TkrTowerAlign;
  using CalibData::TkrTowerAlignCalib;
  using facilities::Util;
  using CLHEP::Hep3Vector;

  std::vector<DOMElement* > towerElts;
  
  Dom::getDescendantsByTagName(docElt, "tower", towerElts);

  // Arg is max tower id
  TkrTowerAlignCalib* pCalib = new TkrTowerAlignCalib(MAX_TOWER_ID);
  refpObject = pCalib;
  
  setBaseInfo(pCalib);

  for (unsigned iElt = 0; iElt < towerElts.size(); iElt++) {
    DOMElement* elt = towerElts[iElt];
    Hep3Vector* pDisp=0;
    Hep3Vector* pRot=0;
    unsigned id;
    StatusCode ret = parseAlignElt(elt, id, pDisp, pRot);
    if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
    if (id > MAX_TOWER_ID)       return StatusCode::FAILURE;
    pCalib->putAlign((unsigned) id, *pDisp, *pRot);
  }
  return StatusCode::SUCCESS;
}

// Now same for TkrInternalAlign

class XmlTkrInternalAlignCnv : public XmlBaseCnv {
  /// Friend needed for instantiation
  friend class CnvFactory<XmlTkrInternalAlignCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();
protected:

  XmlTkrInternalAlignCnv(ISvcLocator* svcs);

  virtual ~XmlTkrInternalAlignCnv() {}       // most likely nothing to do 

  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject);
};

//static CnvFactory<XmlTkrInternalAlignCnv> s_internalFactory;
//const  ICnvFactory& XmlTkrInternalAlignCnvFactory = s_internalFactory;
DECLARE_CONVERTER_FACTORY(XmlTkrInternalAlignCnv);


XmlTkrInternalAlignCnv::XmlTkrInternalAlignCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_TKR_InternalAlign) { 
}

const CLID& XmlTkrInternalAlignCnv::objType() const {
  return CLID_Calib_TKR_InternalAlign;
}

const CLID& XmlTkrInternalAlignCnv::classID() {
  return CLID_Calib_TKR_InternalAlign;
}


StatusCode XmlTkrInternalAlignCnv::i_createObj(const DOMElement* docElt, 
                                               DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::TkrInternalAlign;
  using CalibData::TkrInternalAlignCalib;
  using facilities::Util;
  using CLHEP::Hep3Vector;

  std::vector<DOMElement* > towerElts;

  Dom::getDescendantsByTagName(docElt, "towerInternal", towerElts);

  TkrInternalAlignCalib* pCalib = new TkrInternalAlignCalib();
  refpObject = pCalib;

  setBaseInfo(pCalib);
  
  for (unsigned iElt = 0; iElt < towerElts.size(); iElt++) {
    DOMElement* elt = towerElts[iElt];
    unsigned towerId;
    // Only info associated with tower elt. is its id
    try {
      int intAtt = Dom::getIntAttribute(elt, "id");
      if (intAtt < 0) return StatusCode::FAILURE;
      towerId = intAtt;
    }    catch (std::exception ex) { return StatusCode::FAILURE; }
    if (!pCalib->checkLimits(&towerId)) return StatusCode::FAILURE;

    std::vector<DOMElement* > trayElts;
    Dom::getChildrenByTagName(elt, "tray", trayElts);
    for (unsigned iTrayElt = 0; iTrayElt < trayElts.size(); iTrayElt++) {
      DOMElement* trayElt = trayElts[iTrayElt];
      Hep3Vector* pDisp=0;
      Hep3Vector* pRot=0;
      unsigned trayId;
      StatusCode ret = parseAlignElt(trayElt, trayId, pDisp, pRot);
      if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
      if (!pCalib->checkLimits(&towerId, &trayId)) return StatusCode::FAILURE;
      ret = pCalib->putAlign(pCalib->makeTrayKey(towerId, trayId), 
                             *pDisp, *pRot);
      if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;

      std::vector<DOMElement* > faceElts;
      Dom::getChildrenByTagName(trayElt, "face", faceElts);
      for (unsigned iFaceElt = 0; iFaceElt < faceElts.size(); iFaceElt++) {
        DOMElement* faceElt = faceElts[iFaceElt];
        pDisp = pRot = 0;
        unsigned faceId;
        ret = parseAlignElt(faceElt, faceId, pDisp, pRot);
        if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
        if (!pCalib->checkLimits(&towerId, &trayId, &faceId)) 
          return StatusCode::FAILURE;
        ret = pCalib->putAlign(pCalib->makeFaceKey(towerId, trayId, faceId), 
                               *pDisp, *pRot);
        if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
    
        std::vector<DOMElement* > ladderElts;
        Dom::getChildrenByTagName(faceElt, "ladder", ladderElts);
        for (unsigned iLadderElt = 0; iLadderElt < ladderElts.size(); 
             iLadderElt++) {
          DOMElement* ladderElt = ladderElts[iLadderElt];
          pDisp = pRot = 0;
          unsigned ladderId;
          ret = parseAlignElt(ladderElt, ladderId, pDisp, pRot);
          if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
          if (!pCalib->checkLimits(&towerId, &trayId, &faceId, &ladderId)) 
            return StatusCode::FAILURE;
          ret = pCalib->putAlign(pCalib->makeLadderKey(towerId, trayId, faceId,
                                                       ladderId),*pDisp,*pRot);
          if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;

          std::vector<DOMElement* > waferElts;
          Dom::getChildrenByTagName(ladderElt, "wafer", waferElts);
          for (unsigned iWaferElt = 0; iWaferElt < waferElts.size(); 
               iWaferElt++) {
            pDisp = pRot = 0;
            unsigned waferId;
            ret = parseAlignElt(waferElts[iWaferElt], waferId, pDisp, pRot);
            if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
            if (!pCalib->checkLimits(&towerId, &trayId, &faceId, &ladderId,
                                     &waferId)) return StatusCode::FAILURE;
            ret = pCalib->putAlign(pCalib->makeWaferKey(towerId, trayId,faceId,
                                                        ladderId, waferId), 
                                   *pDisp, *pRot);
            if (ret != StatusCode::SUCCESS) return StatusCode::FAILURE;
          }
        }
      }
    }
  }
  return StatusCode::SUCCESS;
}

