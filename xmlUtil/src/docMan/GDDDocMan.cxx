// $Header$
#include "xmlUtil/docMan/GDDDocMan.h"
#include <xercesc/dom/DOM_Element.hpp>
#include "xmlUtil/Constants.h"
#include "xml/Dom.h"

namespace xmlUtil {
  GDDDocMan* GDDDocMan::m_me = 0;

  GDDDocMan* GDDDocMan::getPointer() {
    if (m_me == 0) {
      m_me = new GDDDocMan();
    }
    return m_me;
  }

  bool GDDDocMan::parse(const std::string& filename, bool saveNotes,
                        const std::string& docType) {
    using namespace xml;

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

  void GDDDocMan::ConstsClient::handleChild(DOM_Node node) {
    if (node == DOM_Node()) { //clean up call
      m_doc = 0;
      if (m_sub != 0) {
        delete m_sub;
        m_sub = 0;
      }
      return;
    }
    if (node.getNodeType() == DOM_Node::DOCUMENT_NODE) { 
      DOM_Document& doc = static_cast<DOM_Document&>(node);

      // Save CVSid and DTDversion attributes
      CVSid = xml::Dom::getAttribute(doc.getDocumentElement(), "CVSid");
      DTDversion = xml::Dom::getAttribute(doc.getDocumentElement(), 
                                          "DTDversion");
      // Evaluate constants
      Constants constants(doc);
      constants.normalizePrimary();
      constants.evalConstants();
      constants.pruneConstants(m_saveNotes);

      m_doc = doc;
      m_sub = new  Substitute(m_doc);
    }  else {  // Child element of root; substitute
      DOM_Element& elt = static_cast<DOM_Element&>(node);
      if (elt.getAttributeNode("substituted") != DOM_Node()) {
        m_sub->execute(elt);
      }
    }

    return;
  }
}      
      
