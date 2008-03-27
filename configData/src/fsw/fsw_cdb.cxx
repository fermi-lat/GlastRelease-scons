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
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

// other headers
#include "configData/base/ConfigXmlUtil.h"
#include "./fsw_cdb_print.h"
#include "./fsw_datum.h"

namespace configData {

  fsw_cdb* fsw_cdb::read_cdb( fsw_enums::CDM_TYPE cdmType, const char* name, int option ) {
    fsw_cdb* retVal(0);
    switch ( cdmType ) {
    case fsw_enums::GFC:
      retVal = new fsw_cdb_inst<GFC_DB_schema>(name,option);
      break;
    case fsw_enums::DFC:
      retVal = new fsw_cdb_inst<DFC_DB_schema>(name,option);
      break;
    case fsw_enums::MFC:
      retVal = new fsw_cdb_inst<MFC_DB_schema>(name,option);
      break;
    case fsw_enums::HFC:
      retVal = new fsw_cdb_inst<HFC_DB_schema>(name,option);
      break;
    default:
      break;
    }
    return retVal;
  }

  
  fsw_enums::CDM_TYPE fsw_cdb::cdm_type( const std::string& cdb_type_name ) {
    static std::map<std::string,fsw_enums::CDM_TYPE> typeNameMap;
    if ( typeNameMap.size() == 0 ) {
      typeNameMap["GFC"] = fsw_enums::GFC;
      typeNameMap["HFC"] = fsw_enums::HFC;
      typeNameMap["DFC"] = fsw_enums::DFC;
      typeNameMap["MFC"] = fsw_enums::MFC;
    }
    static std::map<std::string,fsw_enums::CDM_TYPE>::const_iterator itrFind = typeNameMap.find(cdb_type_name);
    return itrFind != typeNameMap.end() ? itrFind->second : fsw_enums::INVALID;
  }

  fsw_cdb::fsw_cdb(const char* cdm_name)
    :m_name(cdm_name){
  }
  
  fsw_cdb::~fsw_cdb() {
  }  
  
  // print to a stream
  // return 0 for success, error flag otherwise
  int fsw_cdb::printToStream( std::ostream& os ) const {
    os << std::endl;
    os << std::endl;
    os << "File: " << get_name() << std::endl;
    os << "CDM identification" << std::endl;
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::hex << "Key" << std::right
       << ':' << m_key << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::hex << "Schema" << std::right
       << ':' << m_schemaID << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << std::dec << "Version" << std::right
       << ':' << m_versionID << std::endl;   
    os << std::setw(2) << ' ' 
       << std::setw(33) << std::left << "Instance" << std::right
       << ':' << m_instanceID << std::endl;   
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    io->print(0,os);
    return 0;
  }
  
  // write to an XML file
  // return 0 for success, error flag otherwise
  int fsw_cdb::writeToXmlFile( const char* fileName ) const {    

    configData::XmlUtil::init();
    
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    DOMElement* doc = configData::XmlUtil::makeDocument("FSW");
    if ( doc == 0 ) return 2;
    DOMElement* check = io->writeToXml(*doc);
    if ( check == 0 ) return 3;
    configData::XmlUtil::writeIt(*doc,fileName);
    return 0;
  }
  
  // write from an XML file
  // return 0 for success, error flag otherwise    
  int fsw_cdb::readFromXmlFile( const char* fileName ) {
    fsw_datum* io = get_io_handler();
    if ( io == 0 ) return 1;
    //DOMElement* doc = configData::XmlUtil::makeDocument("FSW");
    DOMElement* doc(0);
    if ( doc == 0 ) return 2;
    bool ok = io->readFromXml(*doc);
    if ( ! ok ) return 3;
    return 0;
  }


  template <typename SCHEMA>
  fsw_cdb_inst<SCHEMA>::fsw_cdb_inst(const char* cdm_name, int option)
    :fsw_cdb(cdm_name),
     m_schema(0){
    m_schema = load( option );
    if ( m_schema == 0 ) {
      std::cerr << "Failed to load CDM from " << cdm_name << std::endl;
    } else {
      std::cout << "Loaded " << cdm_name << std::endl;
    }
  }
    
  template <typename SCHEMA>
  fsw_cdb_inst<SCHEMA>::~fsw_cdb_inst() {
    delete m_schema;
  }
  
  template <typename SCHEMA>
  SCHEMA* fsw_cdb_inst<SCHEMA>::load( int option ) {    
    const void* schema_ptr = fsw_get_data( get_name().c_str(), option,
					   &m_key, &m_schemaID, &m_versionID, &m_instanceID );
    m_schema = const_cast<SCHEMA*>(static_cast<const SCHEMA*>(schema_ptr));
    return m_schema;
  }
  
  template <typename SCHEMA>
  fsw_datum* fsw_cdb_inst<SCHEMA>::get_io_handler() const {
    if ( m_schema == 0 ) return 0;
    fsw_datum* retVal = new fsw_datum_inst<SCHEMA>( "Config", *m_schema );
    return retVal;
  }

}
