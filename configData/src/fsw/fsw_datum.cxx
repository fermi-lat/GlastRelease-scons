#ifndef fsw_datum_CXX
#define fsw_datum_CXX
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Base class for converting FSW headers to XML
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
#include "./fsw_datum.h"

// c++/stl headers
#include <iostream>
#include <iomanip>

// Local headers
#include "configData/base/ConfigXmlUtil.h"
#include "./fsw_cdb.h"



namespace configData {

  // write bookkeeping information to an XML file
  DOMElement* fsw_datum::writeXmlHeader( DOMElement& parent, const char* cdm_name, const fmx_id& fmxId ) {
    DOMElement* node = configData::XmlUtil::makeChildNode(parent,"FMX");
    if ( node == 0 ) return 0;
    DOMElement* cnode = configData::XmlUtil::makeChildNodeWithContent(*node,"FMX_FILE",cdm_name);
    if ( cnode == 0 ) return 0;   
    char tmp[16];
    sprintf(tmp,"0x%08x",fmxId.get_key());
    cnode = configData::XmlUtil::makeChildNodeWithContent(*node,"FMX_KEY",tmp);
    if ( cnode == 0 ) return 0;
    sprintf(tmp,"0x%04x",fmxId.get_schemaID());    
    cnode = configData::XmlUtil::makeChildNodeWithContent(*node,"FMX_SCEHMA",tmp);
    if ( cnode == 0 ) return 0;
    cnode = configData::XmlUtil::makeChildNodeWithContent(*node,"FMX_VERSION",fmxId.get_versionID());
    if ( cnode == 0 ) return 0;
    cnode = configData::XmlUtil::makeChildNodeWithContent(*node,"FMX_INSTANCE",fmxId.get_instanceID());
    if ( cnode == 0 ) return 0;    
    return node;
  }

  // C'tor for for fsw_datum gives a name to any item of data
  fsw_datum::fsw_datum(const char* name) 
    :m_name(name){
  }
  
  // D'tor is trivial
  fsw_datum::~fsw_datum() {
  }
  
  
  // C'tor for fsw_datum_mask, gives a name to a 32 bit mask
  fsw_datum_mask::fsw_datum_mask(const char* name, unsigned int& val)
    :fsw_datum(name),
     m_mask(val){
  }

  // D'tor is trivial
  fsw_datum_mask::~fsw_datum_mask() {
  }


  // Print the mask to a stream
  void fsw_datum_mask::print( int depth, std::ostream& os ) const {
    
    // print the name
    os << std::setw(depth) << ' '  
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ": ";    
    // print 32 bits
    unsigned int test(0x1);
    for ( int i(0); i < 32; i++ ) {
      os << ( m_mask & test ? '.' : 'X' );
      test *= 2;
    }
    os << std::endl;
  }

  // Write the mask to a stream
  DOMElement* fsw_datum_mask::writeToXml( DOMElement& parent ) const {
    // use format "0x%08x", ie. zero-filled 8 digit hex preceded by "0x"
    char tmp[12];
    sprintf(tmp,"0x%08x",m_mask);
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), tmp);
  }
  
  // Read the mask from a stream
  bool fsw_datum_mask::readFromXml( DOMElement& node ) {
    return true;
  }  

}

#endif
