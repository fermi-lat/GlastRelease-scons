// $Header

#include "xmlUtil/id/DictField.h"
#include "xml/Dom.h"

namespace xmlUtil {
  DictField::DictField(DOM_Element elt) {
    DOM_Element child = xml::Dom::getFirstChildElement(elt);

    if (child != DOM_Element())
      m_constraints = new DictConstraints(child);
    else m_constraints = 0;

    m_name = std::string(xml::Dom::transToChar(elt.getAttribute("name")));
  }

  bool DictField::allowed(const unsigned value) {
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
