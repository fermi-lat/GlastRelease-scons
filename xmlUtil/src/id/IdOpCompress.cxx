// $Header$

#include "xmlUtil/id/IdOpCompress.h"
#include "xml/Dom.h"

namespace xmlUtil {
  IdOpCompress::IdOpCompress(DOM_Element cmprss) {
    from = xml::Dom::getAttribute(cmprss, "from");
    to = xml::Dom::getAttribute(cmprss, "to");
  }

  NamedId *IdOpCompress::convert(const NamedId& inputId) {

    int fromIx = inputId.hasField(from);
    int toIx = inputId.hasField(to);

    // For now if input id doesn't have "from" and "to" fields
    // in the proper order, return null pointer.  Another
    // possible response in this case would be to return a copy
    // of the input id.

    if (toIx < fromIx) return 0;
    if (fromIx < 0) return 0;

    NamedId *newId = new NamedId(inputId, toIx + 1);

    NamedId::Fields *newFields = newId->m_fields;
    NamedId::Fields *inFields = inputId.m_fields;

    for (unsigned int ix = fromIx; ix < inFields->size(); ix++ ) {
      NamedId::IdField * aField = new NamedId::IdField( *((*inFields)[ix]) );
      newFields->push_back(aField);
    }
    return newId;
  }

}
