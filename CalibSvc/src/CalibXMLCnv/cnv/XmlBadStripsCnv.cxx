// $Header$

#include <string>
#include "XmlBadStripsCnv.h"
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

// #include "CalibData/Tkr/BadStrips.h"  already in XmlBadStripsCnv.h
#include "CalibData/CalibTime.h"
#include "xml/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

static CnvFactory<XmlBadStripsCnv> s_factory;
const  ICnvFactory& XmlBadStripsCnvFactory = s_factory;

XmlBadStripsCnv::XmlBadStripsCnv( ISvcLocator* svc) :
  XmlBaseCnv(svc, CLID_Calib_TKR_BadChan) {
}

const CLID& XmlBadStripsCnv::objType() const {
  return CLID_Calib_TKR_BadChan;
}

const CLID& XmlBadStripsCnv::classID() {
  return CLID_Calib_TKR_BadChan;
}

namespace {


  /** Internal helper.  Convert string list of integers to entries in vector
    @param s      Input string of one or more pos. integers separated by spaces
    @param strips Vector of short unsigned, output

  */
  void strToNum(const std::string& s, 
                CalibData::StripCol* v){
  
    std::string::const_iterator it = s.begin();
    
    // Maybe add something to be sure all we've got are digits
    // and blanks??
    
    // Skip over leading blanks, if any
    while((it != s.end()) && (*it == ' ')) it++;
    
    // save contiguous digits in tempStr
    while ((it != s.end()) && (*it >= '0') && (*it <= '9'))    {
      std::string tempStr;    
      tempStr += *it;
      it++;  
      
      while((it != s.end()) && (*it >= '0') && (*it <= '9')){
        tempStr += *it;
        it++;
      }
      
      // No more contiguous digits; skip over blanks
      while((it != s.end()) && (*it == ' ')) it++;
      
      // Save the converted integer
      v->push_back(atoi(tempStr.c_str()));
    }
  }
}


// Create our specific object, given
StatusCode XmlBadStripsCnv::i_createObj(const DOM_Element& element,
                                        DataObject*& refpObject)
{
  using xml::Dom;
  using CalibData::BadStrips;

  // Interesting information in file is the bad type (hot or dead)
  // and per-tower information.
  std::string typeAtt;
  typeAtt = Dom::getAttribute(element, "badType");

  BadStrips::eBadType bType;
  if (typeAtt.compare("hot") == 0) bType = BadStrips::HOT;
  else if (typeAtt.compare("dead") == 0) bType = BadStrips::DEAD;
  else {
    // complain?
    return StatusCode::FAILURE;
  }
  BadStrips* pBad = 
    new BadStrips(bType, *m_vstart, *m_vend, m_serNo);
  refpObject = pBad;

  // Find child tower elements 
  
  // const 
  StatusCode sc = StatusCode::FAILURE;
  DOM_Element child = Dom::findFirstChildByName(element, "tower");
  while (child != DOM_Element()) {
    sc = processTower(child, pBad);

    // if bad return do something

    child = Dom::getSiblingElement(child);
  }

  return sc;
}

StatusCode XmlBadStripsCnv::processTower(const DOM_Element& towerElt,
                                         CalibData::BadStrips *pBad) {
  using xml::Dom;

  std::string attValue;
  attValue = Dom::getAttribute(towerElt, "row");
  unsigned row = (unsigned) atoi(attValue.c_str());

  attValue = Dom::getAttribute(towerElt, "col");
  unsigned col = (unsigned) atoi(attValue.c_str());

  attValue = Dom::getAttribute(towerElt, "allBad");
  bool allBad = (attValue.compare("true") == 0);

  int howBad = 0;
  if (allBad) {
    attValue = Dom::getAttribute(towerElt, "howBad");
    howBad = atoi(attValue.c_str());

  } 

  StatusCode sc = 
    pBad->addBadTower(allBad, howBad, row, col);

  if (allBad) return sc;

  // Process each uniplane element
  DOM_Element uniElt = Dom::getFirstChildElement(towerElt);

  while (uniElt != DOM_Element()) {
    sc = processUni(uniElt, row, col, pBad);
    // if bad status, complain and return
    uniElt = Dom::getSiblingElement(uniElt);
  }
  return sc;
                                  
}

StatusCode XmlBadStripsCnv::processUni(const DOM_Element& uniElt, 
                                       unsigned row, unsigned col,
                                       CalibData::BadStrips* pBad)
{
  using CalibData::BadStrips;
  using xml::Dom;

  CalibData::StripCol strips;

  std::string attValue;
  
  attValue = Dom::getAttribute(uniElt, "allBad");
  bool allBad = (attValue.compare("true") == 0);

  attValue = Dom::getAttribute(uniElt, "howBad");
  int howBad = atoi(attValue.c_str());

  attValue = Dom::getAttribute(uniElt, "tray");
  unsigned int tray = (unsigned int)atoi(attValue.c_str());

  attValue = Dom::getAttribute(uniElt, "which");
  bool top;
  if (attValue.compare("top") == 0) top = true;
  else if (attValue.compare("bot") == 0) top = false;
  // if anything else happens, xml file is bad

  if (!allBad) {   // process strip lists, spans
    addStrips(uniElt, &strips);
  }

  // finally,
  pBad->addBadPlane(row, col, tray, top, howBad, strips);
 
  return StatusCode::SUCCESS;
}

StatusCode XmlBadStripsCnv::addStrips(const DOM_Element& uniElt,
                                      CalibData::StripCol* strips) {
  using xml::Dom;

  // Children of uniElt are
  //    at most one stripList, followed by
  //    arbitrary number of span elements
  DOM_Element child = Dom::findFirstChildByName(uniElt, "stripList");
  if (child != DOM_Element()) {
    std::string xmlList = Dom::getAttribute(child, "strips");
    strToNum(xmlList, strips);
    child = Dom::getSiblingElement(child);
  }

  while (child != DOM_Element()) {
    std::string firstStr = Dom::getAttribute(child, "first");
    unsigned short first = (unsigned short) atoi(firstStr.c_str());
    std::string lastStr = Dom::getAttribute(child, "last");
    unsigned short last = (unsigned short) atoi(lastStr.c_str());
      
    if (last >= first) {
      // Might as well reserve memory all at once
      strips->reserve(strips->size() + last + 1 - first);  
      for (unsigned short int i = first; i <= last; i++) {
        strips->push_back(i);
      }
    }
    child = Dom::getSiblingElement(child);
  }
  return StatusCode::SUCCESS;
}

