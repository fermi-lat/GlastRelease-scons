// $Header$

#include "xmlUtil/id/IdConversion.h"
#include "xmlUtil/id/IdOperation.h"
#include "xml/Dom.h"
#include "dom/DOMString.hpp"
#include "xmlUtil/id/IdOpTruncate.h"
#include "xmlUtil/id/IdOpDisappear.h"
#include "xmlUtil/id/IdOpCompress.h"

namespace xmlUtil {
  IdConversion::IdConversion() : m_path(0), m_condition(0) {
    m_op = new IdOperation();
  }

  IdConversion::IdConversion(DOM_Element conversion) {

    // Get first child; invoke private function to build path
    DOM_Element child = xml::Dom::getFirstChildElement(conversion);
    makePath(child);

    // Get next child; save field name in condition component
    child = xml::Dom::getSiblingElement(child);
    m_condition = new std::string(xml::Dom::getAttribute(child, "name"));

    // Get next child;  build new op component.
    child = xml::Dom::getSiblingElement(child);
    buildOp(child);
  }

  bool IdConversion::inDomain(const NamedId& inputId) {
    return inputId.hasSubpath(*m_path);
  }

  bool IdConversion::satisfies(const NamedId& inputId) {
    if (m_condition == 0) return true;

    return (inputId.hasField(*m_condition) >= 0);
  }
    
  NamedId * IdConversion::convert(const NamedId& inputId) {
    if (!inDomain(inputId)) return 0;
    return internalConvert(inputId);
  }

  NamedId * IdConversion::internalConvert(const NamedId& inputId) {
    if (satisfies(inputId)) {  // let the operation do its thing
      return m_op->convert(inputId);
    }
    else { // clone input and return
      return new NamedId(inputId);
    }
  }

  void IdConversion::makePath(const DOM_Element& pathElt) {
    
    // "path" consists of a list of fields.  Fields have 
    // a required attribute "name".  Save its value.
    m_path = new NameSeq;
    DOM_Element child = xml::Dom::getFirstChildElement(pathElt);

    while (child != DOM_Element()) {
      m_path->push_back(new std::string(xml::Dom::getAttribute(child, "name")));
      child = xml::Dom::getSiblingElement(child);
    }

  }

  // Could maybe do something more elegant than a switch, but it would
  // take quite a bit of machinery and would only be worthwhile if
  // the set of ops was expected to change often.
  void IdConversion::buildOp(const DOM_Element& opElt) {
    DOMString opType = opElt.getTagName();

    if (opType.equals("truncate"))
    {
      m_op = new IdOpTruncate(opElt);
    }
    else if (opType.equals("disappear")) {
      m_op = new IdOpDisappear(opElt);
    }
    else if (opType.equals(DOMString("compress"))) {
      m_op = new IdOpCompress(opElt);
    }
    else { // default to identity operation, implemented by base class
      m_op = new IdOperation(opElt);
    }
  }

  // return true if our path is subpath of given conversion's path
  bool IdConversion::subpathOf(const IdConversion& other) const {
    unsigned int ourLen = m_path->size();
    if (ourLen > other.m_path->size() ) return false;

    for (unsigned int ix = 0; ix < ourLen; ix++) {
      if ((*(other.m_path))[ix]->compare((*(*m_path)[ix])) ) return false;
    }
    return true;
  }

  ostream& operator<<(ostream& s, const IdConversion& convers) {
    s << (*(convers.m_op)) << std::endl << " Path: " << (*(convers.m_path)) << 
      " Condition: hasField " << (*(convers.m_condition));
    return s;
  }

}
