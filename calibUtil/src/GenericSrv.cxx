//$ Header: $
/// Module implements methods for clients to get generic services.

#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>

#include "calibUtil/GenericSrv.h"

namespace calibUtil {
 
  /// Constructor that fills in values of generic data from XML file
  GenericSrv::GenericSrv(const DOM_Element& docElt) : m_sample(0) {
    
    DOM_Element child  = xml::Dom::findFirstChildByName(docElt,"generic");
    m_instName  = xml::Dom::getAttribute(child,"instrument");
    m_timestamp = xml::Dom::getAttribute(child,"timestamp");
    m_calType   = xml::Dom::getAttribute(child,"calType");
    m_fmtVer    = xml::Dom::getAttribute(child,"fmtVersion");

  }
  
}// end of namespace calibUtil






