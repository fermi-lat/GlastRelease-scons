// $Header$

#ifndef XMLUTIL_IDOPCOMPRESS_H
#define XMLUTIL_IDOPCOMPRESS_H

#include <string>
#include "id/IdOperation.h"


namespace xmlUtil {
  //!
  class IdOpCompress : public IdOperation {
  public:
    IdOpCompress(const DOM_Element) {};
    ~IdOpCompress() {};

    virtual *NamedId convert(const NamedId& inputId);
  private:
    std:string from;
    std:string to;

  }    

}
#endif
