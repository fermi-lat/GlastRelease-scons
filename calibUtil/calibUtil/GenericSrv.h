// $Header$
 
#ifndef CALIBUTIL_GENERICSRV_H
#define CALIBUTIL_GENERICSRV_H

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"

#include <string>
#include <iostream>
#include <vector>

namespace calibUtil {

  class InputSample;

  class GenericSrv {
  public: 

    /// Constructor that fills in values of generic data from XML file
    GenericSrv(const DOM_Element& docElt);

    /// Constructor to be used when creating a new calibration data set
    GenericSrv(std::string inst, std::string timestamp, 
               std::string calType, InputSample* samp) :
      m_instName(inst), m_timestamp(timestamp), m_calType(calType),
      m_fmtVer(0), m_sample(samp) {}

    /// Get instrument name
    std::string getInst(){
      return m_instName;
    }
    
    /// Get timestamp
    std::string getTimestamp(){
      return m_timestamp;
    }
    
    /// Get calibration type
    std::string getCalType(){
      return m_calType;
    }
    
    /// Get format Version
    std::string getFmtVer(){
      return m_fmtVer;
    }
    
  private:
    
    std::string m_instName;
    std::string m_timestamp;
    std::string m_calType;
    std::string m_fmtVer;

    // Corresponds to inputSample element in badStrips.dtd (except
    // that the class doesn't actually exist yet.
    InputSample *m_sample;

  };

}// end of namespace calibUtil

#endif



