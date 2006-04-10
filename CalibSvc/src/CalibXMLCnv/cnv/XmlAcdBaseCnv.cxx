// $Header$
#include "GaudiKernel/MsgStream.h"
#include "XmlAcdBaseCnv.h"
#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMNode.hpp>
#include "facilities/Util.h"

// #include "idents/AcdId.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

XmlAcdBaseCnv::XmlAcdBaseCnv(ISvcLocator* svc, const CLID& clid) :
  XmlBaseCnv(svc, clid),  m_nDet(10000), m_nPmt(10000), m_id(0) {}

StatusCode XmlAcdBaseCnv::readAcdDimension(const DOMElement* docElt, 
                                        unsigned& nDet, unsigned& nFace,
                                        unsigned& nRow, unsigned& nCol, 
                                        unsigned& nNA, unsigned& nPmt) {
  using xmlBase::Dom;

  MsgStream log(msgSvc(), "XmlBaseCnv" );
  DOMElement* dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == 0) return StatusCode::FAILURE;

  try {
    nDet = Dom::getIntAttribute(dimElt, "nTile");
    nFace = Dom::getIntAttribute(dimElt, "nFace");
    nRow = Dom::getIntAttribute(dimElt, "nRow");
    nCol = Dom::getIntAttribute(dimElt, "nCol");
    nPmt = Dom::getIntAttribute(dimElt, "nPmt");
    nNA = Dom::getIntAttribute(dimElt, "nNA");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::readAcdDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }
  return StatusCode::SUCCESS;
}

DOMElement* XmlAcdBaseCnv::findFirstPmt(const DOMElement* docElt) {
  using xmlBase::Dom;

  // note that so-called <tile> may be ribbon or NA
  DOMElement* tileElt = Dom::findFirstChildByName(docElt, "tile");
  // If no <tile> elements, this file is useless
  if (tileElt == 0) return tileElt;
  DOMElement* pmtElt = Dom::getFirstChildElement(tileElt);
  if (pmtElt == 0) return pmtElt;

  try {
    std::string idString = Dom::getAttribute(tileElt, "tileId");
    m_id = idents::AcdId(idString);
    m_nPmt = Dom::getIntAttribute(pmtElt, "iPmt");
  }
  catch (xmlBase::DomException ex1) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::findFirstPmt" << std::endl;
    std::cerr << ex1.getMsg() << std::endl;
    throw ex1;
  }
  catch (facilities::WrongType ex2) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::findFirstPmt" << std::endl;
    std::cerr << ex2.what() << std::endl;
    throw ex2;
  }

  return pmtElt;
}

/// Still another one to navigate XML file and find next set of pmt data
DOMElement* XmlAcdBaseCnv::findNextPmt(const DOMElement* pmtElt) {
  using xmlBase::Dom;


  DOMElement* elt = Dom::getSiblingElement(pmtElt);

  if (elt != 0) {
    m_nPmt = Dom::getIntAttribute(elt, "iPmt");
    return elt;
  }

  // Done with this pmt; look for sibling
  DOMNode* node = pmtElt->getParentNode();
  DOMElement* tileElt = static_cast<DOMElement* >(node);   // current tile
  tileElt = Dom::getSiblingElement(tileElt);          // next tile

  if (tileElt != 0) {
    try {
      std::string idString = Dom::getAttribute(tileElt, "tileId");
      m_id = idents::AcdId(idString);

      elt = Dom::getFirstChildElement(tileElt);
      m_nPmt = Dom::getIntAttribute(elt, "iPmt");
    }
    catch (xmlBase::DomException ex1) {
      std::cerr << "From CalibSvc::XmlAcdBaseCnv::findNextPmt" << std::endl;
      std::cerr << ex1.getMsg() << std::endl;
      throw ex1;
    }
    catch (facilities::WrongType ex2) {
      std::cerr << "From CalibSvc::XmlAcdBaseCnv::findNextPmt" << std::endl;
      std::cerr << ex2.what() << std::endl;
      throw ex2;
    }
  }
  return elt;
}
