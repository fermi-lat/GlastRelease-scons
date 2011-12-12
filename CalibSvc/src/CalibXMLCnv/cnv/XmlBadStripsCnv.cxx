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
#include "xmlBase/Dom.h"

// For channel status bit definitions
#include "calibUtil/ChannelStatusDef.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlBadStripsCnv> s_factory;
//const  ICnvFactory& XmlBadStripsCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlBadStripsCnv);

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
StatusCode XmlBadStripsCnv::i_createObj(const DOMElement* element,
                                        DataObject*& refpObject)
{
  using xmlBase::Dom;
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
  StatusCode sc = StatusCode::SUCCESS;   // it's OK to have no bad towers
  DOMElement* child = Dom::findFirstChildByName(element, "tower");
  while (child != 0) {
    sc = processTower(child, pBad);

    // if bad return do something

    child = Dom::getSiblingElement(child);
  }

  return sc;
}

StatusCode XmlBadStripsCnv::processTower(const DOMElement* towerElt,
                                         CalibData::BadStrips *pBad) {
  using xmlBase::Dom;

  //  std::string attValue;
  unsigned row, col;
  try {
    row = Dom::getIntAttribute(towerElt, "row");
    col = Dom::getIntAttribute(towerElt, "col");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlBadStripsCnv::processTower" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  bool allBad = 0;
  int howBad = 0;
  std::string attValue;

  attValue = Dom::getAttribute(towerElt, "nOnbdCalib");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdCalib;
  }
  attValue = Dom::getAttribute(towerElt, "nOnbdTrig");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdTrig;
  }
  attValue = Dom::getAttribute(towerElt, "nOnbdData");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdData;
  }
  allBad = (howBad != 0);

  StatusCode sc = 
    pBad->addBadTower(allBad, howBad, row, col);

  if (allBad) return sc;

  // Process each uniplane element
  DOMElement* uniElt = Dom::getFirstChildElement(towerElt);

  while (uniElt != 0) {
    sc = processUni(uniElt, row, col, pBad);
    // if bad status, complain and return
    uniElt = Dom::getSiblingElement(uniElt);
  }
  return sc;
                                  
}

StatusCode XmlBadStripsCnv::processUni(const DOMElement* uniElt, 
                                       unsigned row, unsigned col,
                                       CalibData::BadStrips* pBad)
{
  using CalibData::BadStrips;
  using xmlBase::Dom;

  CalibData::StripCol strips;

  std::string attValue;
  
  attValue = Dom::getAttribute(uniElt, "allBad");
  bool allBad = (attValue.compare("true") == 0);

  int howBad = 0;

  attValue = Dom::getAttribute(uniElt, "nOnbdCalib");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdCalib;
  }
  attValue = Dom::getAttribute(uniElt, "nOnbdTrig");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdTrig;
  }
  attValue = Dom::getAttribute(uniElt, "nOnbdData");
  if (attValue.compare("true") == 0) {
    howBad |= vCALIBUTIL_nOnbdData;
  }

  unsigned int tray;
  try {    
    //  attValue = Dom::getAttribute(uniElt, "tray");
    //  unsigned int tray = (unsigned int)atoi(attValue.c_str());
    tray = Dom::getIntAttribute(uniElt, "tray");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlBadStripsCnv::processUni" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  attValue = Dom::getAttribute(uniElt, "which");
  bool top;
  if (attValue.compare("top") == 0) top = true;
  else if (attValue.compare("bot") == 0) top = false;
  // if anything else happens, xml file is bad

  if (!allBad) {   // process strip lists, spans
    addStrips(uniElt, &strips);
  }

  // finally,
  pBad->addBadPlane(row, col, tray, top, howBad, allBad, strips);
 
  return StatusCode::SUCCESS;
}

StatusCode XmlBadStripsCnv::addStrips(const DOMElement* uniElt,
                                      CalibData::StripCol* strips) {
  using xmlBase::Dom;

  // Children of uniElt are an arbitrary collection of stripList
  // and stripSpan elements
  DOMElement* childElt = Dom::getFirstChildElement(uniElt);

  while (childElt != 0 ) {
    // must be list or span
    //    if ((childElt.getTagName()).equals("stripList") ) {
    if (Dom::checkTagName(childElt, "stripList")) {
      std::string xmlList = Dom::getAttribute(childElt, "strips");
      strToNum(xmlList, strips);
    }
    else if (Dom::checkTagName(childElt, "stripSpan")) {
        //    else if ((childElt.getTagName()).equals("stripSpan") ) {
      unsigned short first, last;
      try {
        first = Dom::getIntAttribute(childElt, "first");
        last = Dom::getIntAttribute(childElt, "last");
      }
      catch  (xmlBase::DomException ex) {
        std::cerr << "From CalibSvc::XmlBadStripsCnv::addStrips" << std::endl;
        std::cerr << ex.getMsg() << std::endl;
        throw ex;
      }

      if (last >= first) {
        // Might as well reserve memory all at once
        strips->reserve(strips->size() + last + 1 - first);  
        for (unsigned short int i = first; i <= last; i++) {
          strips->push_back(i);
        }
      }
    }
    childElt = Dom::getSiblingElement(childElt);
  }

  return StatusCode::SUCCESS;
}

