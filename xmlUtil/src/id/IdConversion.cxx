// $Header$

#include "xmlUtil/id/IdConversion.h"
#include "xmlUtil/id/IdOperation.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/id/IdOpTruncate.h"
#include "xmlUtil/id/IdOpDisappear.h"
#include "xmlUtil/id/IdOpCompress.h"

namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE
  IdConversion::IdConversion() : m_path(0), m_condition(0) {
    m_op = new IdOperation();
  }

  IdConversion::IdConversion(const DOMElement* conversion) {

    // Get first child; invoke private function to build path
    DOMElement* child = xmlBase::Dom::getFirstChildElement(conversion);
    makePath(child);

    // Get next child; save field name in condition component
    child = xmlBase::Dom::getSiblingElement(child);
    m_condition = new std::string(xmlBase::Dom::getAttribute(child, "name"));

    // Get next child;  build new op component.
    child = xmlBase::Dom::getSiblingElement(child);
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

  void IdConversion::makePath(const DOMElement* pathElt) {
    
    // "path" consists of a list of fields.  Fields have 
    // a required attribute "name".  Save its value.
    m_path = new NameSeq;
    DOMElement* child = xmlBase::Dom::getFirstChildElement(pathElt);

    while (child != 0) {
      m_path->push_back(new std::string(xmlBase::Dom::getAttribute(child, "name")));
      child = xmlBase::Dom::getSiblingElement(child);
    }

  }

  // Could maybe do something more elegant than a switch, but it would
  // take quite a bit of machinery and would only be worthwhile if
  // the set of ops was expected to change often.
  void IdConversion::buildOp(const DOMElement* opElt) {

    std::string opType = xmlBase::Dom::getTagName(opElt);

    if (opType == std::string("truncate"))
    {
      m_op = new IdOpTruncate(opElt);
    }
    else if (opType == std::string("disappear")) {
      m_op = new IdOpDisappear(opElt);
    }
    else if (opType == std::string("compress")) {
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

  std::ostream& operator<<(std::ostream& s, const IdConversion& convers) {
    s << (*(convers.m_op)) << std::endl << " Path: " << (*(convers.m_path)) << 
      " Condition: hasField " << (*(convers.m_condition));
    return s;
  }

}
