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
  GenericSrv::GenericSrv(DOM_Element docElt){
    
    DOM_Element child  = xml::Dom::findFirstChildByName(docElt,"generic");
    instName  = xml::Dom::getAttribute(child,"instrument");
    timestamp = xml::Dom::getAttribute(child,"timestamp");
    calType   = xml::Dom::getAttribute(child,"calType");
    fmtVer    = xml::Dom::getAttribute(child,"fmtVersion");

  }
  
}// end of namespace calibUtil






