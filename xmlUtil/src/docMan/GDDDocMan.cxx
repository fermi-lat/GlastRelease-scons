// $Header$
#include "xmlUtil/docMan/GDDDocMan.h"
#include "dom/DOM_Element.hpp"
#include "xmlUtil/Constants.h"

namespace xmlUtil {
  GDDDocMan* GDDDocMan::m_me = 0;

  GDDDocMan* GDDDocMan::getPointer() {
    if (m_me == 0) {
      m_me = new GDDDocMan();
    }
    return m_me;
  }

  bool GDDDocMan::parse(const std::string& filename, bool saveNotes) {
    using namespace xml;

    m_constsClient->m_saveNotes = saveNotes;
    return DocMan::parse(filename);
  }

  GDDDocMan::GDDDocMan() {
    m_constsClient = new ConstsClient(this);
    regMeFirst(m_constsClient);
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
      // Evaluate constants
      DOM_Document& doc = static_cast<DOM_Document&>(node);
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
      
