// $Header$

#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {
  NamedId::NamedId(const int len) {
    /* For now -- until I understand how to change
       addField to make use of pre-allocation, always
       allocate 0-length vector.
    */
    m_fields = new Fields(0);
  }

  NamedId::~NamedId() {
    delete m_fields;
  }

  NamedId::NamedId(const NamedId& toCopy) {
    m_fields = new Fields(*(toCopy.m_fields));
  }

  void NamedId::addField(const std::string name, const unsigned value) {
    IdField* field = new IdField;
    field->name = name;
    field->value = value;
    m_fields->push_back(field);
  }

  bool NamedId::hasPath(NameSeq& path) {
    int pathLen = path.size();
    if (pathLen > this->size()) return false;

    // Else check one by one
    NamedId::iterator idIt = this->begin();
    NameSeq::iterator pathIt = path.begin();

    for (int ix = 0; ix < pathLen; ix++) {
      if (path[ix].compare((*m_fields)[ix].name) ) return false;
    }
    return true;
  }

  bool NamedId::hasField(const std::string fname) {
    for (int ix = 0; ix < size(); ix++) {
      if (!(fname.compare((*m_fields)[ix].name) ) ) return true;
    }
    return false;
  }
}

