// $Header$

#ifndef XMLUTIL_IDCONVERTERLESSTHAN_H
#define XMLUTIL_IDCONVERTERLESSTHAN_H

#include <functional>
#include "xmlUtil/id/IdConverter.h"

#include "xmlUtil/id/IdConversion.h"

namespace xmlUtil {
  /*! This class describes a function object satisfying the
      definition of Strict Weak Ordering (like < except that
        a <= b and b <= a does not imply a == b)
      so that it may be used by the STL \b sort algorithm
  */
  class IdConverterLessThan : 
    public std::binary_function<IdConverter::ConversionIt, 
    IdConverter::ConversionIt, bool> 
  {
  public:
    bool operator()(IdConversion* x, 
                    IdConversion* y) {
      return ((x->path)->size() < ((y->path)->size() ) );
    }
  };
}

#endif
