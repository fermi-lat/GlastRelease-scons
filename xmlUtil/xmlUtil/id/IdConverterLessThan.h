// $Header$

#ifndef XMLUTIL_IDCONVERTERLESSTHAN_H
#define XMLUTIL_IDCONVERTERLESSTHAN_H

#include "xmlUtil/id/IdConverter.h"

#include "xmlUtil/id/IdConversion.h"

namespace xmlUtil {
  /*! This class describes a function object satisfying the
      definition of Strict Weak Ordering (like < except that
        a <= b and b <= a does not imply a == b)
      so that it may be used by the STL \b sort algorithm
  */
  struct IdConverterLessThan : 
    public binary_function<IdConverter::ConversionIt, 
    IdConverter::ConversionIt, bool> 
  {
    bool operator()(IdConversion* x, 
                    IdConversion* y) {
      return ((x->path)->size() < ((y->path)->size() ) );
    }
  };
}

#endif
