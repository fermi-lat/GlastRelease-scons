// $Header$
#include "GaudiKernel/MsgStream.h"
#include "XmlCalBaseCnv.h"
#include "xmlBase/Dom.h"
#include "idents/CalXtalId.h"

// Local utilities to interpret attributes 
namespace {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

  unsigned findRangeAtt(const DOMElement* rangeElt) {
    using xmlBase::Dom;
    using idents::CalXtalId;

    std::string att = Dom::getAttribute(rangeElt, "range");
    if (att.size() == 0) {    // Try diode next
      att = Dom::getAttribute(rangeElt, "diode");
      if (att.size() == 0)      return 0;      // no sort of range

      // else process value of diode attribute;
      //   a little bit of sleight of hand here
      if (att.compare(std::string("LE")) == 0) return CalXtalId::LEX8;
      if (att.compare(std::string("HE")) == 0) return CalXtalId::HEX8;
      return 0;
    }
    if (att.compare(std::string("LEX8")) == 0) return CalXtalId::LEX8;
    if (att.compare(std::string("LEX1")) == 0) return CalXtalId::LEX1;
    if (att.compare(std::string("HEX8")) == 0) return CalXtalId::HEX8;
    if (att.compare(std::string("HEX1")) == 0) return CalXtalId::HEX1; 
    // anything else is illegal.  Should be caught by parser, but
    // maybe should also throw exception here.
    return 0;
  }

  unsigned findFace(const DOMElement* faceElt) {
    using xmlBase::Dom;
    using idents::CalXtalId;

    std::string att = Dom::getAttribute(faceElt, "end");
    if (att.compare(std::string("NEG")) == 0) return CalXtalId::NEG;
    if (att.compare(std::string("POS")) == 0) return CalXtalId::POS;
    return 0;     // in case "end" is not applicable, this is the answer
  }
}

// redundant since this is already in XmlBaseCnv.h
using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

XmlCalBaseCnv::XmlCalBaseCnv(ISvcLocator* svc, const CLID& clid) :
  XmlBaseCnv(svc, clid), m_nRow(10000), m_nCol(10000), 
  m_nLayer(10000), m_nXtal(10000),
  m_nFace(10000), m_nRange(10000) {}

StatusCode XmlCalBaseCnv::readDimension(const DOMElement* docElt, 
                                        unsigned& nRow, unsigned& nCol, 
                                        unsigned& nLayer,
                                        unsigned& nXtal, unsigned& nFace,
                                        unsigned& nRange,
                                        unsigned* nDacCol,
                                        unsigned* nXpos) {
  using xmlBase::Dom;

  MsgStream log(msgSvc(), "XmlCalBaseCnv" );
  DOMElement* dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == 0) return StatusCode::FAILURE;

  try {
    nRow = Dom::getIntAttribute(dimElt, "nRow");
    nCol = Dom::getIntAttribute(dimElt, "nCol");
    nLayer = Dom::getIntAttribute(dimElt, "nLayer");
    nXtal = Dom::getIntAttribute(dimElt, "nXtal");
    nFace = Dom::getIntAttribute(dimElt, "nFace");
    nRange = Dom::getIntAttribute(dimElt, "nRange");
    if (nDacCol) *nDacCol = Dom::getIntAttribute(dimElt, "nDacCol");
    if (nXpos) *nXpos = Dom::getIntAttribute(dimElt, "nXpos");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlCalBaseCnv::readDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  // Can't check dimensions if <dimension> description isn't meant to be
  // exact.  No "exact" attribute ==> exact = true
  std::string exactAttr = Dom::getAttribute(dimElt, "exact");
  if (exactAttr == std::string("false")) return StatusCode::SUCCESS;

  unsigned expected = nRow * nCol;
  // Make some consistency checks.  # tower elements should be nRow*nCol
  // # layer elements should be nRow*nCol*nLayer, etc.
  std::vector<DOMElement*> nlist;

  Dom::getDescendantsByTagName(docElt, "tower", nlist);
  
  if (nlist.size() != expected) {
    log << MSG::ERROR << "Expected tower dimension <> actual  " << endreq;
    // put out a message and...
    return StatusCode::FAILURE;
  }
  expected *= nLayer;
  Dom::getDescendantsByTagName(docElt, "layer", nlist);
  if (nlist.size() != expected) {
    log << MSG::ERROR << "Expected layer dimension <> actual  " << endreq;
    // put out a message and...
    return StatusCode::FAILURE;
  }

  expected *= nXtal;
  Dom::getDescendantsByTagName(docElt, "xtal", nlist);
  if (nlist.size() != expected) {
    log << MSG::ERROR << "Expected xtal dimension <> actual  " << endreq;
    // put out a message and...
    return StatusCode::FAILURE;
  }
  
  expected *= nFace;
  Dom::getDescendantsByTagName(docElt, "face", nlist);
  if (nlist.size() != expected) {
    log << MSG::ERROR << "Expected face dimension <> actual  " << endreq;
    // put out a message and...
    return StatusCode::FAILURE;
  }

  expected *= nRange;
  // Not as easy to check for the right number of range elements here
  // since they have different names, depending on what type of calibration
  // this is.  However they all have the same name, so just find first
  // child of a face element, then count all similarly-named elements

  DOMElement* child = Dom::getFirstChildElement(nlist[0]);
  std::string tagName = Dom::getTagName(child);
  Dom::getDescendantsByTagName(docElt, tagName, nlist);
  //  nlist = docElt.getElementsByTagName(child.getTagName());
  if (nlist.size() != expected) {
    log << MSG::ERROR << "Expected calib type elements " << expected 
        << " <> actual  " << nlist.size() << endreq;
    // put out a message and...
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

DOMElement* XmlCalBaseCnv::findFirstRange(const DOMElement* docElt) {
  using xmlBase::Dom;
  using idents::CalXtalId;

  DOMElement* elt = Dom::findFirstChildByName(docElt, "tower");
  if (elt == 0) return elt;

  //  std::string att = Dom::getAttribute(elt, "nRow");
  //  m_nRow = atoi(att.c_str());
  try {
    m_nRow = Dom::getIntAttribute(elt, "iRow");
    m_nCol = Dom::getIntAttribute(elt, "iCol");


    // All child elements of a tower are layer elements
    elt = Dom::getFirstChildElement(elt);
    m_nLayer = Dom::getIntAttribute(elt, "iLayer");

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);
    m_nXtal = Dom::getIntAttribute(elt, "iXtal");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlCalBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  // All child elements of xtal are face elements
  elt = Dom::getFirstChildElement(elt); 
  m_nFace = findFace(elt);
  
  elt = Dom::getFirstChildElement(elt);
  m_nRange = findRangeAtt(elt);
  return elt;
}

/// Still another one to navigate XML file and find next set of range data
DOMElement* XmlCalBaseCnv::findNextRange(const DOMElement* rangeElt) {
  using xmlBase::Dom;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

  DOMElement* elt = Dom::getSiblingElement(rangeElt);
  if (elt != 0) {
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this xtal face; look for sibling
  DOMNode* node = rangeElt->getParentNode();
  elt = static_cast<DOMElement* >(node);   // current xtal face
  elt = Dom::getSiblingElement(elt);          // next xtal face

  if (elt != 0) {
    m_nFace = findFace(elt);

    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this xtal
  node = node->getParentNode();  // current xtal
  elt = static_cast<DOMElement* >(node);
  elt = Dom::getSiblingElement(elt);         // next xtal

  if (elt != 0) {
    try {
      m_nXtal = Dom::getIntAttribute(elt, "iXtal");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    //    std::string att = Dom::getAttribute(elt, "iXtal");
    //    m_nXtal = atoi(att.c_str());

    // All child elements of xtal are face elements
    elt = Dom::getFirstChildElement(elt);
    m_nFace = findFace(elt);

    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this layer
  node = node->getParentNode();  // current layer
  elt = static_cast<DOMElement* >(node);
  elt = Dom::getSiblingElement(elt);         // next layer

  if (elt != 0) {
    try {
      m_nLayer = Dom::getIntAttribute(elt, "iLayer");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }
    

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);

    try {
      m_nXtal = Dom::getIntAttribute(elt, "iXtal");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlCalBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    // All child elements of xtal are face elements
    elt = Dom::getFirstChildElement(elt);
    m_nFace = findFace(elt);
    
    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this tower
  node = node->getParentNode();  // current tower
  elt = static_cast<DOMElement* >(node);
  elt = Dom::getSiblingElement(elt);         // next tower

  if (elt == 0) return elt;

  // otherwise we've got a new tower; go through the whole
  // rigamarole


  try {
    m_nRow = Dom::getIntAttribute(elt, "iRow");
    m_nCol = Dom::getIntAttribute(elt, "iCol");


    // All child elements of a tower are layer elements
    elt = Dom::getFirstChildElement(elt);
    m_nLayer = Dom::getIntAttribute(elt, "iLayer");

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);
    m_nXtal = Dom::getIntAttribute(elt, "iXtal");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlCalBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }


  // All child elements of xtal are face elements
  elt = Dom::getFirstChildElement(elt);
  m_nFace = findFace(elt);

  elt = Dom::getFirstChildElement(elt);
  m_nRange = findRangeAtt(elt);

  return elt;
}
