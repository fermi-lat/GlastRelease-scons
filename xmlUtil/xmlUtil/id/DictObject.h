// $Header$
#ifndef XMLUTIL_DICTOBJECT_H
#define XMLUTIL_DICTOBJECT_H

namespace xmlUtil {
  //! Pure virtual base class to support visitors to an id dictionary;
  //! all classes associated with an id dictionary which may be
  //! visited inherit from this class.
  //! Currently, this is all id dictionary classes except
  //! DictConstraints

  class DictVisitor;

  class DictObject {
  public:

    virtual bool accept(DictVisitor* vis) = 0;
    virtual ~DictObject() {}
  };

}
#endif
