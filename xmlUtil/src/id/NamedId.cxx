// $Header$

#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {
  //  typedef std::vector <IdField* > Fields;

  //  typedef Fields::iterator  FieldIt;

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

  //! Just copy first n fields to a new NamedId
  NamedId::NamedId(const NamedId& toCopy, unsigned int n) {
    unsigned int ix;
    m_fields = new Fields(0);

    for (ix = 0; ix < n; ix++) {
      (*m_fields)[ix] = new IdField(*(*(toCopy.m_fields))[ix]);
    }

  }

  void NamedId::addField(std::string name, unsigned value) {
    IdField* field = new IdField;
    field->name = name;
    field->value = value;
    m_fields->push_back(field);
  }

  void NamedId::addField(const IdField& newField) {
    IdField* field = new IdField(newField);
    m_fields->push_back(field);
  }

  bool NamedId::hasPath(const NameSeq& path) const {
    unsigned int pathLen = path.size();
    if (pathLen > m_fields->size()) return false;

    // Else check one by one
    //    FieldIt idIt = m_fields->begin();
    //    NameSeq::iterator pathIt = path.begin();

    for (unsigned int ix = 0; ix < pathLen; ix++) {
      if (path[ix]->compare(((*m_fields)[ix])->name) ) return false;
    }
    return true;
  }

  int NamedId::hasField(std::string fname) const {
    for (unsigned int ix = 0; ix < m_fields->size(); ix++) {
      if (!(fname.compare(((*m_fields)[ix])->name) ) ) return ix;
    }
    return -1;
  }

}

