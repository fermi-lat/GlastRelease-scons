//$ Header: $
/// Module implements methods for clients to get generic services.

#include "xml/Dom.h"

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>

#include "calibUtil/GenericSrv.h"

namespace calibUtil {
 
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  /// Constructor that fills in values of generic data from XML file
  GenericSrv::GenericSrv(const DOMElement* docElt) : m_sample(0) {
    
    DOMElement* child  = xml::Dom::findFirstChildByName(docElt,"generic");
    m_instName  = xml::Dom::getAttribute(child,"instrument");
    m_timestamp = xml::Dom::getAttribute(child,"timestamp");
    m_calType   = xml::Dom::getAttribute(child,"calType");
    m_fmtVer    = xml::Dom::getAttribute(child,"fmtVersion");

  }
  
}// end of namespace calibUtil






