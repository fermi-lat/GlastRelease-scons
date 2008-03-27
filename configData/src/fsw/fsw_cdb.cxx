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
#include "./fsw_cdb.h"

// c++/stl headers
#include <iostream>
#include <iomanip>
#include <fstream>

// other headers
#include "configData/base/ConfigXmlUtil.h"
#include "./fsw_datum.h"


namespace configData {
  
  // C'tor for fsw_cdb
  fsw_cdb::fsw_cdb(const char* cdm_name, 
		   unsigned int key, short schemaID, 
		   short versionID, short instanceID)
    :m_fmxId(key,schemaID,versionID,instanceID),
     m_name(cdm_name){     
  }

  /// Dtor for fsw_cdb
  fsw_cdb::~fsw_cdb() {
  }  
  

  // print to a stream
  int fsw_cdb::printToStream( std::ostream& os ) const {
    // Get some space
    os << std::endl << std::endl;
    // Print bookkeeping stuff
    os << "File: " << get_name() << std::endl;
    os << "CDM identification" << std::endl;
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::hex << "Key" << std::right
       << ':' << m_fmxId.get_key() << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::hex << "Schema" << std::right
       << ':' << m_fmxId.get_schemaID() << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::dec << "Version" << std::right
       << ':' << m_fmxId.get_versionID() << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << "Instance" << std::right
       << ':' << m_fmxId.get_instanceID() << std::endl;   
    // Get the templated object that know about the structure
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    // print the CDM using the handler
    io->print(0,os);
    return 0;
  }
  

  // write to an XML file
  int fsw_cdb::writeToXmlFile( const char* fileName ) const {            
    // Get the templated object that know about the structure
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    // Make sure that Xerces is ready
    configData::XmlUtil::init();
    // Make the TOP level document
    DOMElement* doc = configData::XmlUtil::makeDocument("FSW");
    if ( doc == 0 ) return 2;    
    DOMElement* node = configData::fsw_datum::writeXmlHeader(*doc,m_name.c_str(),m_fmxId);
    if ( node == 0 ) return 2;
    // convert the CDM to XML using the handler
    node = io->writeToXml(*doc);
    if ( node == 0 ) return 2;
    // Write the XML
    if ( ! XmlUtil::writeIt(*doc,fileName) ) return 3;
    return 0;
  }
  
  // read from an XML file
  int fsw_cdb::readFromXmlFile( const char* fileName ) {
    // Get the templated object that know about the structure     
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    //DOMElement* doc = configData::XmlUtil::makeDocument("FSW");
    DOMElement* doc(0);
    if ( doc == 0 ) return 2;
    bool ok = io->readFromXml(*doc);
    if ( ! ok ) return 3;
    return 0;
  }

}

