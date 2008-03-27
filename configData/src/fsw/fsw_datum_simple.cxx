#ifndef fsw_datum_simple_CXX
#define fsw_datum_simple_CXX
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Partial specializations for handling IO of simple types in CDMs
//
// Environment:
//      Software developed for GLAST.
//
// Author List:
//      Eric Charles
//
// Copyright Information:
//      Copyright (C) 2008      Stanford Linear Accelerator Center
//
//---------------------------------------------------------------------------


// base class headers
#include "./fsw_datum_inst.h"

// c++/stl headers
#include <iostream>
#include <iomanip>

// Local headers
#include "configData/base/ConfigXmlUtil.h"

// Other headers
#include "facilities/Util.h"

namespace configData {

  // Print a pointer to a stream
  template <>
  void fsw_datum_inst<void*>::specializedPrint( int depth, std::ostream& os) const {
    // just print the pointer value
    os << std::setw(depth) << ' '  
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ": " << m_datum << std::endl;
  }

  // Print an unsigned int to a stream
  template <>
  void fsw_datum_inst<unsigned int>::specializedPrint( int depth, std::ostream& os) const {
    // just print the value
    os << std::setw(depth) << ' ' 
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ": " << m_datum << std::endl;
  }


  // Write a pointer xml
  template <>
  DOMElement* fsw_datum_inst<void*>::specializedWriteToXml( DOMElement& parent ) const {
    // looks like <var_name>0x_ptrVal</varName>
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), m_datum);
  }

  // Write an unsigned int to xml
  template <>
  DOMElement* fsw_datum_inst<unsigned int>::specializedWriteToXml( DOMElement& parent ) const {
    // looks like <var_name>value</varName>
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), m_datum );
  }

  // Read a pointer from XML
  template <>
  bool fsw_datum_inst<void*>::specializedReadFromXml( DOMElement& thisNode ) {
    // This is a no-op
    return true;
  }
  
  // Read an unsigned int from XML
  template <>
  bool fsw_datum_inst<unsigned int>::specializedReadFromXml( DOMElement& thisNode ) {
    std::string txtVal;
    XmlUtil::getText(thisNode,txtVal);
    if ( txtVal.empty() ) return false;
    try { m_datum = facilities::Util::stringToUnsigned(txtVal); }
    catch ( std::exception& ) { return false; }
    return true;
  }

}

#endif
