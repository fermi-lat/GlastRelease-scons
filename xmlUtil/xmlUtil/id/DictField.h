// $Header$

#ifndef XMLUTIL_DICTFIELD_H
#define XMLUTIL_DICTFIELD_H

#include <string>

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/DictConstraints.h"

namespace xmlUtil {
  //! Represent a field -- its name plus constraints on its values, if any --
  //! which may be used as part of a NamedId

  class DictField {
  public: 
    DictField(DOM_Element elt);
    ~DictField(){ if (m_constraints) delete m_constraints;};

    const std::string& getName() const {return m_name;}

    //! Does the specified value satisfy the constraints, if any?
    bool allowed(const unsigned value);

  private:
    friend class IdDict;  // may also need friend IdConverter
    friend class DictNode;

    // Build a field which has value constraints
    DictField(const std::string& name, const DictConstraints& constraints)
      :   m_name(name) 
    {m_constraints = new DictConstraints(constraints);}

    // Build a field with no value constraints
    DictField(const std::string& name) : m_constraints(0), m_name(name) {}

    DictConstraints* getConstraints() const {return m_constraints;}

    DictConstraints *m_constraints;
    std::string     m_name;
    
  };       // end DictField class
}          // end xmlUtil namespace

#endif
