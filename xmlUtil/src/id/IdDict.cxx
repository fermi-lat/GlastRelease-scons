// $Header$

#include "dom/DOMString.hpp"
#include "dom/DOM_NodeList.hpp"
#include "xml/Dom.h"
#include "xmlUtil/id/IdDict.h"
#include "xmlUtil/id/DictFieldMan.h"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictNode.h"

namespace xmlUtil {
  IdDict::IdDict(DOM_Element elt) {
    // Check that element has the right tag name: idDict
    // Caller probably will have done this already
    if (!(elt.getTagName().equals("idDict"))){ ;} // big trouble

    // Check number of children.  This is an upper bound
    // on number of fields we have to store
    int size = (elt.getChildNodes()).getLength();
    m_fieldMan = new DictFieldMan(size);
                
    DOM_Element fieldElt = xml::Dom::getFirstChildElement(elt);

    // Register all the fields
    while ((fieldElt.getNodeName()).equals("field")) {
      //      DOM_Element fieldElt = static_cast<DOM_Element&>(fieldNode);
      DictField* field = new DictField(fieldElt);
      m_fieldMan->signup(field);
      fieldElt = xml::Dom::getSiblingElement(fieldElt);
    }
    // Finally make the hierarchy of constraints on Identifiers
    m_root = new DictNode(fieldElt, 0, m_fieldMan);
  }
  
  bool IdDict::addChild(DictNode* parent, DictNode* newNode) {
    /* Started to write something to add field if not already
       registered, but whoever has created newNode should
       already have done this */
    return parent->addChild(newNode);
  }
}
