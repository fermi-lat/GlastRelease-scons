// $Header

#include "xmlUtil/id/DictField.h"
#include "xml/Dom.h"

namespace xmlUtil {
  DictField::DictField(DOM_Element elt) {
    DOM_Element child = xml::Dom::getFirstChildElement(elt);

    if (child != DOM_Element())
      m_constraints = new DictConstraints(child);
    else m_constraints = 0;
    std::string name(xml::Dom::transToChar(child.getAttribute("name")));
    m_name = name;
  }

  bool DictField::allowed(const unsigned value) {
    if (m_constraints) return m_constraints->allowed(value);
    else return true;
  }

}
