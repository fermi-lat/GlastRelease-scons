// $Header $

#include "xmlUtil/id/IdOperation.h"

namespace xmlUtil {

  // The only reason for putting this here is that it generates
  // a warning message.  Better one warning message than many.
  IdOperation::IdOperation(DOM_Element) {}

  NamedId *IdOperation::convert(const NamedId& inputId)      {
    return new NamedId(inputId);
  }

  ostream& operator<<(ostream& s, const IdOperation& op) {
    s << op.myOp();
    return s;
  }
}
