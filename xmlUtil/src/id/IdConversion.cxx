// $Header$

#include "xmlUtil/id/IdConversion.h"
#include "xmlUtil/id/IdOperation.h"

namespace xmlUtil {
  IdConversion::IdConversion() : path(0), condition(0) {
    op = new IdOperation();
  }

  IdConversion::IdConversion(const DOM_Element) {

    /*   ******** TO DO ************ */

  }

  bool IdConversion::inDomain(NamedId& inputId) {
    return inputId.hasPath(*path);
  }

  bool IdConversion::satisfies(NamedId& inputId) {
    if (condition == 0) return true;

    return inputId.hasField(*condition);
  }
    
  *NamedId IdConversion::convert(NamedId& inputId) {
    if (!inDomain(inputId)) return 0;
    return internalConvert(inputId);
  }

  *NamedId IdConversion::internalConvert(NamedId& inputId) {
    if (satisfies(inputId)) {  // let the operation do its thing
      return op->convert(inputId);
    }
    else { // clone input and return
      return new NamedId(inputId);
    }
  }





}
