// $Header$

#include "xmlUtil/id/IdOperation.h"

namespace xmlUtil {

  // The only reason for putting this here is that it generates
  // a warning message.  Better one warning message than many.
  IdOperation::IdOperation(const XERCES_CPP_NAMESPACE_QUALIFIER DOMElement*) {}

  NamedId *IdOperation::convert(const NamedId& inputId)      {
    return new NamedId(inputId);
  }

  std::ostream& operator<<(std::ostream& s, const IdOperation& op) {
    s << op.myOp();
    return s;
  }
}
