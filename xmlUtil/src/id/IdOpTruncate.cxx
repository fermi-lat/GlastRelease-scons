// $Header$

#include "xmlUtil/id/IdOpTruncate.h"
#include "xmlBase/Dom.h"

namespace xmlUtil {

  IdOpTruncate::IdOpTruncate(const XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* trunc)
  {
    // Get values of attributes "start" and "beyond", save
    start = xmlBase::Dom::getAttribute(trunc, "start");
    std::string beyondVal = xmlBase::Dom::getAttribute(trunc, "beyond");
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
