// $Header$

#include "xmlUtil/id/IdOpTruncate.h"
#include "xml/Dom.h"

namespace xmlUtil {

  IdOpTruncate::IdOpTruncate(DOM_Element trunc) {
    // Get values of attributes "start" and "beyond", save
    start = xml::Dom::getAttribute(trunc, "start");
    std::string beyondVal = xml::Dom::getAttribute(trunc, "beyond");
    beyond = (beyondVal.compare(std::string("true")) == 0);
  }

  NamedId * IdOpTruncate::convert(const NamedId& inputId) {
    // Check that inputId has a "start" field.  
    int nToCopy;
    if ((nToCopy = inputId.hasField(start)) < 0) return 0;

    if (beyond) ++nToCopy;
    NamedId *newId = new NamedId(inputId, nToCopy);
    return newId;
  }
}
