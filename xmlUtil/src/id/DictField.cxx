// $Header

#include "xmlUtil/id/DictField.h"
#include "xml/Dom.h"

namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE
  DictField::DictField(DOMElement* elt) {
    DOMElement* child = xml::Dom::getFirstChildElement(elt);

    if (child != 0) m_constraints = new DictConstraints(child);
    else m_constraints = 0;

    m_name = xml::Dom::getAttribute(elt, "name");
      //std::string(xml::Dom::transToChar(elt->getAttribute("name")));
  }

  bool DictField::allowed(unsigned value) const {
    if (m_constraints) return m_constraints->allowed(value);
    else return true;
  }

  DictField::DictField(const DictField& toCopy) : DictObject(toCopy) {
    deepCopy(toCopy);
  }

  DictField& DictField::operator=(const DictField& d) {
    if (this != &d) {
      if (m_constraints) delete m_constraints;
      deepCopy(d);
    }
    return *this;
  }

  void DictField::deepCopy(const DictField& toCopy) {
    m_name = toCopy.m_name;
    if (toCopy.m_constraints) {
      m_constraints = new DictConstraints(*(toCopy.m_constraints));
    }
    else m_constraints = 0;
  }

}
