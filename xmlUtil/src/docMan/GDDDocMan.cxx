// $Header$
#include "xmlUtil/docMan/GDDDocMan.h"
#include <xercesc/dom/DOMElement.hpp>
#include "xmlUtil/Constants.h"
#include "xmlBase/Dom.h"

namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE
  GDDDocMan* GDDDocMan::m_me = 0;

  GDDDocMan* GDDDocMan::getPointer() {
    if (m_me == 0) {
      m_me = new GDDDocMan();
    }
    return m_me;
  }

  bool GDDDocMan::parse(const std::string& filename, bool saveNotes,
                        const std::string& docType) {
    using namespace xmlBase;

    m_constsClient->m_saveNotes = saveNotes;
    return DocMan::parse(filename, docType);
  }

  GDDDocMan::GDDDocMan() {
    m_constsClient = new ConstsClient(this);
    regMeFirst(m_constsClient);
  }

  const std::string& GDDDocMan::getCVSid() const {
    return m_constsClient->getCVSid();
  }
  const std::string& GDDDocMan::getDTDversion() const {
    return m_constsClient->getDTDversion();
  }

  void GDDDocMan::ConstsClient::handleChild(XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* node) {
    if (node == 0) { //clean up call
      m_doc = 0;
      if (m_sub != 0) {
        delete m_sub;
        m_sub = 0;
      }
      return;
    }
    if (node->getNodeType() == DOMNode::DOCUMENT_NODE) { 
      // Evaluate constants
      //      DOMDocument* doc = static_cast<DOMDocument *>(node);
      DOMDocument* doc = static_cast<DOMDocument *>(node);
      CVSid = xmlBase::Dom::getAttribute(doc->getDocumentElement(), "CVSid");
      DTDversion = xmlBase::Dom::getAttribute(doc->getDocumentElement(), 
                                          "DTDversion");
      Constants constants(doc);
      constants.normalizePrimary();
      constants.evalConstants();
      constants.pruneConstants(m_saveNotes);

      m_doc = doc;
      m_sub = new  Substitute(m_doc);
    }  else {  // Child element of root; substitute
      DOMElement* elt = static_cast<DOMElement *>(node);
      //      if (elt->getAttributeNode("substituted") != 0) {
      if (xmlBase::Dom::hasAttribute(elt, "substituted")) {
        m_sub->execute(elt);
      }
    }

    return;
  }
}      
      
