// $Header$

#ifndef XMLUTIL_IDCONVERSION_H
#define XMLUTIL_IDCONVERSION_H

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {
  class DictField;
  class IdOperation;

  //! An /b IdConversion object specifies a set of NamedIds and a
  //! collection of operations it will perform on any of them to
  //! produce a new NamedId.
  class IdConversion {
  public:
    // Default constructor produces a conversion which acts on any NamedId 
    // and does nothing to it; that is, returns a copy.
    IdConverion() {};
    IdConversion(const DOM_Element); 

    ~IdConversion();

    *NamedId convert(NamedId& inputId);
    // {return new NamedId(*inputId);};

  private:
    typedef vector<DictField*> Path;

    Path *path;

    typedef vector<IdOperation*> Operations;

    Operations *ops;
  };  // end of class IdConversion
}   // end of namespace
