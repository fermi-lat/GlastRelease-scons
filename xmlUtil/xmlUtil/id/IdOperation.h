// $Header$

#ifndef XMLUTIL_IDOPERATION_H
#define XMLUTIL_IDOPERATION_H

#include "dom/DOM_Element.hpp"

namespace xmlUtil {
  class NameId;


  //! An /b IdOperation object specifies a transformation which may
  //! be applied to a NamedId object by invoking IdOperation::convert
  /*! Derived classes will in general need some private data to
    keep track of the operation to be performed.
  */
  class IdOperation {
  public:

    //! Will be overridden by derived classes to actually do something
    IdOperation(const DOM_Element) {}; 

    ~IdOperation();

    //! Convert input NamedId in place.
    /*! Return true if /b inputId has whatever fields might be
        needed for the operation.  This base class just leaves the
        the input alone and always returns true. */
    bool convert(NamedId* inputId) {return true;};

    // Default constructor produces an operation which acts on any NamedId 
    // and does nothing to it.  
    // Not yet clear whether this can be made private
    IdOperation() {};
  };
}
