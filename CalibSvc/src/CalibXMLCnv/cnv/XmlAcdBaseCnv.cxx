// $Header$
#include "GaudiKernel/MsgStream.h"
#include "XmlAcdBaseCnv.h"
#include "xmlBase/Dom.h"
#include <xercesc/dom/DOMNode.hpp>

// #include "idents/AcdId.h"

using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

XmlAcdBaseCnv::XmlAcdBaseCnv(ISvcLocator* svc, const CLID& clid) :
  XmlBaseCnv(svc, clid), m_nFace(10000), m_nRow(10000), m_nCol(10000), 
  m_nPmt(10000), m_nRange(10000) {}

StatusCode XmlAcdBaseCnv::readAcdDimension(const DOMElement* docElt, 
                                        unsigned& nFace,
                                        unsigned& nRow, unsigned& nCol, 
                                        unsigned& nPmt, unsigned& nRange) {
  using xmlBase::Dom;

  MsgStream log(msgSvc(), "XmlBaseCnv" );
  DOMElement* dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == 0) return StatusCode::FAILURE;

  try {
    nFace = Dom::getIntAttribute(dimElt, "nFace");
    nRow = Dom::getIntAttribute(dimElt, "nRow");
    nCol = Dom::getIntAttribute(dimElt, "nCol");
    nPmt = Dom::getIntAttribute(dimElt, "nPmt");
    nRange = Dom::getIntAttribute(dimElt, "nRange");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::readAcdDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }
  return StatusCode::SUCCESS;
}

DOMElement* XmlAcdBaseCnv::findFirstRange(const DOMElement* docElt) {
  using xmlBase::Dom;
  //  using idents::CalXtalId;

  /*  Rewrite for ACD */
  DOMElement* faceElt = Dom::findFirstChildByName(docElt, "face");
  // If no <face> elements, this file is useless
  if (faceElt == 0) return faceElt;
  DOMElement* rowElt = Dom::getFirstChildElement(faceElt);
  if (rowElt == 0) return rowElt;
  DOMElement* colElt = Dom::getFirstChildElement(rowElt);
  if (colElt == 0) return colElt;
  DOMElement* pmtElt = Dom::getFirstChildElement(colElt);
  if (pmtElt == 0) return pmtElt;

  DOMElement* rangeElt = Dom::getFirstChildElement(pmtElt);
  if (rangeElt == 0) return rangeElt;
  
  try {
    m_nFace = Dom::getIntAttribute(faceElt, "iFace");
    m_nRow = Dom::getIntAttribute(rowElt, "iRow");
    m_nCol = Dom::getIntAttribute(colElt, "iCol");
    m_nPmt = Dom::getIntAttribute(pmtElt, "iPmt");
    m_nRange = Dom::getIntAttribute(rangeElt, "range");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  return rangeElt;
}

/// Still another one to navigate XML file and find next set of range data
DOMElement* XmlAcdBaseCnv::findNextRange(const DOMElement* rangeElt) {
  using xmlBase::Dom;


  DOMElement* elt = Dom::getSiblingElement(rangeElt);

  if (elt != 0) {
    m_nRange = Dom::getIntAttribute(elt, "range");
    return elt;
  }


  // Done with this pmt; look for sibling
  DOMNode* node = rangeElt->getParentNode();
  elt = static_cast<DOMElement* >(node);   // current pmt
  elt = Dom::getSiblingElement(elt);          // next pmt

  if (elt != 0) {
    m_nPmt = Dom::getIntAttribute(elt, "iPmt");

    elt = Dom::getFirstChildElement(elt);
    m_nRange = Dom::getIntAttribute(elt, "range");
    return elt;
  }


  // Done with this <col>
  node = node->getParentNode();  // current <col> element
  elt = static_cast<DOMElement* >(node);
  elt = Dom::getSiblingElement(elt);         // next <col>

  if (elt != 0) {
    try {
      m_nCol = Dom::getIntAttribute(elt, "iCol");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    // All child elements of col are pmt elements
    elt = Dom::getFirstChildElement(elt);
    m_nPmt = Dom::getIntAttribute(elt, "iPmt");

    elt = Dom::getFirstChildElement(elt);
    m_nRange = Dom::getIntAttribute(elt, "range");
    return elt;
  }

  // Done with this row
  node = node->getParentNode();  // current row
  elt = static_cast<DOMElement* >(node);
  elt = Dom::getSiblingElement(elt);         // next row

  if (elt != 0) {              // find first range in row
    try {
      m_nRow = Dom::getIntAttribute(elt, "iRow");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }
    

    // All child elements of row elements are column elements
    elt = Dom::getFirstChildElement(elt);

    try {
      m_nCol = Dom::getIntAttribute(elt, "iCol");
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    // All child elements of col are pmt elements
    elt = Dom::getFirstChildElement(elt);
    m_nPmt = Dom::getIntAttribute(elt, "iPmt");
    
    elt = Dom::getFirstChildElement(elt);
    m_nRange = Dom::getIntAttribute(elt, "range");
    return elt;
  }

  // Done with this face
  node = node->getParentNode();  // current face
  elt = static_cast<DOMElement*>(node);
  elt = Dom::getSiblingElement(elt);         // next face

  if (elt == 0) return elt;

  // otherwise we've got a new face; go through the whole
  // rigamarole


  try {
    m_nFace = Dom::getIntAttribute(elt, "iFace");

    // All child elements of a face are row elements
    elt = Dom::getFirstChildElement(elt);
    m_nRow = Dom::getIntAttribute(elt, "iRow");

    // All child elements of row elements are col elements
    elt = Dom::getFirstChildElement(elt);
    m_nCol = Dom::getIntAttribute(elt, "iCol");


    // All child elements of col are pmt elements
    elt = Dom::getFirstChildElement(elt);
    m_nPmt = Dom::getIntAttribute(elt, "iPmt");
  }
  catch (xmlBase::DomException ex) {
    std::cerr << "From CalibSvc::XmlAcdBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  // All child elements of pmt elements are the thing we actually want!
  elt = Dom::getFirstChildElement(elt);
  m_nRange = Dom::getIntAttribute(elt, "range");

  return elt;
}
