// $Header$
 
#ifndef CALIBUTIL_GENERICSRV_H
#define CALIBUTIL_GENERICSRV_H

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"

#include <string>
#include <iostream>
#include <vector>

namespace calibUtil {

  class GenericSrv {
  public: 

    /// Constructor that fills in values of generic data from XML file
    GenericSrv(DOM_Element docElt){

      DOM_Element child  = xml::Dom::findFirstChildByName(docElt,"generic");
      instName  = xml::Dom::getAttribute(child,"instrument");
      timestamp = xml::Dom::getAttribute(child,"timestamp");
      calType   = xml::Dom::getAttribute(child,"calType");
      fmtVer    = xml::Dom::getAttribute(child,"fmtVersion");

    }

    /// Get instrument name
    std::string getInst(){
      return instName;
    }

    /// Get timestamp
    std::string getTimestamp(){
      return timestamp;
    }

    /// Get calibration type
    std::string getCalType(){
      return calType;
    }
    
    /// Get format Version
    std::string getFmtVer(){
      return fmtVer;
    }

  private:
    
  std::string instName;
  std::string timestamp;
  std::string calType;
  std::string fmtVer;

  };

}/// end of namespace calibUtil

#endif



