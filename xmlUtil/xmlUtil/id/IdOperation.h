// $Header$

#ifndef XMLUTIL_IDOPERATION_H
#define XMLUTIL_IDOPERATION_H

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/NamedId.h"

namespace xmlUtil {


  //! An /b IdOperation object specifies a transformation which may
  //! be applied to a NamedId object by invoking IdOperation::convert
  /*! Logically IdOperation is a nested class of IdConversion, but
    it's too messy to implement that way. 
    Derived classes will in general need some private data to
    keep track of the operation to be performed.
    Concrete derived classes will likely include
         - disappear   (return null)
         - truncate    (eliminate trailing id fields from specified place)
         - compress    (eliminate id fields between specified endpoints)
  */
  class IdOperation {
  public:

    virtual ~IdOperation() {};
    friend class IdConversion;  //! Must be a friend to invoke convert

    // Perhaps need to throw exception if inputId isn't compatible
    // with requested operation??

  protected:
    // Convert input NamedId in place.  Do we really want this?? 
    // Probably not.
    /* Return true if /b inputId has whatever fields might be
        needed for the operation.  This base class just leaves the
        the input alone and always returns true. */
    // virtual bool convert(NamedId* inputId) {return true;};

    //! Return a pointer to the new id, if any;
    //! behavior of base class is just to copy
    virtual NamedId * convert(const NamedId& inputId);
    // Might also want something to convert a simple identifier (i.e.,
    // just a list of integers without accompanying field names).
    // However converters work in terms of fields names, so one needs
    // an Id Dictionary w.r.t which one can generate a NamedId out of
    // the Identifier, then invoke the above version of convert.

    //! Will be overridden by derived classes to actually do something
    IdOperation(DOM_Element);

    // Default constructor produces an operation which acts on any NamedId 
    // and does nothing to it.  
    // Not yet clear whether this can be made private
    IdOperation() {}
  };
}
#endif

