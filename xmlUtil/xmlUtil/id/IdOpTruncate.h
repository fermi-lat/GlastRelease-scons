// $Header$

#ifndef XMLUTIL_IDOPTRUNCATE_H
#define XMLUTIL_IDOPTRUNCATE_H

#include <string>
#include "xmlUtil/id/IdOperation.h"


namespace xmlUtil {
  //! 
  class IdOpTruncate : public IdOperation {
  public:
    IdOpTruncate(DOM_Element trunc);
    ~IdOpTruncate() {}

    virtual NamedId * convert(const NamedId& inputId);
  private:
    std::string start;
    bool       beyond;
  };    
}
#endif
