// $Header$
#include "GaudiKernel/MsgStream.h"
#include "XmlAncBaseCnv.h"
#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMNode.hpp>
#include "facilities/Util.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

XmlAncBaseCnv::XmlAncBaseCnv(ISvcLocator* svc, const CLID& clid) :
  XmlBaseCnv(svc, clid),  m_iMod(10000), m_iLay(10000), m_iChan(10000) {}

StatusCode XmlAncBaseCnv::readAncDimension(const DOMElement* docElt, 
                                           unsigned& nMod, unsigned& nLay,
                                           unsigned& nChan) {
  using xmlBase::Dom;

  MsgStream log(msgSvc(), "XmlBaseCnv" );
  DOMElement* dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == 0) return StatusCode::FAILURE;

  try {
    nMod = Dom::getIntAttribute(dimElt, "nModule");
    nLay = Dom::getIntAttribute(dimElt, "nLayer");
    nChan = Dom::getIntAttribute(dimElt, "nChannel");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlAncBaseCnv::readAncDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }
  m_isTagger = m_hasLayers = (nLay > 0);
  m_isQdc = !m_isTagger;
  if (m_isTagger) m_eltName = "taggerChan";
  if (m_isQdc) {
    m_eltName = "qdcChan";
    m_iLay = 0;
  }
  return StatusCode::SUCCESS;
}

DOMElement* XmlAncBaseCnv::findFirstChan(const DOMElement* docElt) {
  using xmlBase::Dom;

  // note that so-called <tile> may be ribbon or NA
  DOMElement* chanElt = Dom::findFirstChildByName(docElt, m_eltName);

  if (chanElt == 0) return chanElt;

  try {
    m_iMod = Dom::getIntAttribute(chanElt, "iMod");
    m_iChan = Dom::getIntAttribute(chanElt, "iChan");
    if (m_hasLayers) {
      m_iLay = Dom::getIntAttribute(chanElt, "iLay");
    }
  }
  catch (xmlBase::DomException ex1) {
    std::cerr << "From CalibSvc::XmlAncBaseCnv::findFirstChan" << std::endl;
    std::cerr << ex1.getMsg() << std::endl;
    throw ex1;
  }
  catch (facilities::WrongType ex2) {
    std::cerr << "From CalibSvc::XmlAncBaseCnv::findFirstChan" << std::endl;
    std::cerr << ex2.what() << std::endl;
    throw ex2;
  }

  return chanElt;
}

/// Still another one to navigate XML file and find next set of chan data
/// Save attributes indicating which channel it is.
DOMElement* XmlAncBaseCnv::findNextChan(const DOMElement* elt) {
  using xmlBase::Dom;

  DOMElement* chanElt = Dom::getSiblingElement(elt);
  if (chanElt == 0) return chanElt;

  try {
    m_iMod = Dom::getIntAttribute(chanElt, "iMod");
    m_iChan = Dom::getIntAttribute(chanElt, "iChan");
    if (m_hasLayers) {
      m_iLay = Dom::getIntAttribute(chanElt, "iLay");
    }
    return chanElt;
  }
  catch (xmlBase::DomException ex1) {
    std::cerr << "From CalibSvc::XmlAncBaseCnv::findNextChan" << std::endl;
    std::cerr << ex1.getMsg() << std::endl;
    throw ex1;
  }
  catch (facilities::WrongType ex2) {
    std::cerr << "From CalibSvc::XmlAncBaseCnv::findNextChan" << std::endl;
    std::cerr << ex2.what() << std::endl;
    throw ex2;
  }
  return 0;
}
