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
#include "./fsw_datum.h"

// FSW headers
#include "fsw/GFC_DB/GFC_DB_schema.h"
#include "fsw/XFC_DB/HFC_DB_schema.h"
#include "fsw/XFC_DB/MFC_DB_schema.h"
#include "fsw/XFC_DB/DFC_DB_schema.h"



namespace configData {

  fsw_cdb* fsw_cdb::read_cdb( const char* cdm_name, int option ) {

    unsigned int key(0);
    short schemaID(-1);
    short versionID(-1);
    short instanceID(-1);

    void* schema_ptr = const_cast<void*>(fsw_get_data( cdm_name, option,
						       &key, &schemaID, &versionID, &instanceID ));

    fsw_cdb* retVal(0);
    switch ( schemaID ) {
    case fsw_enums::GFC_SCHEMA_ID :
      retVal = new fsw_cdb_inst<GFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<GFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::DFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<DFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<DFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::MFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<MFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<MFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::HFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<HFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<HFC_DB_schema*>(schema_ptr));
      break;
    default:
      break;
    }
    return retVal;
  } 

  fsw_cdb::fsw_cdb(const char* cdm_name, 
		   unsigned int key, short schemaID, 
		   short versionID, short instanceID)
    :m_key(key),
     m_schemaID(schemaID),
     m_versionID(versionID),
     m_instanceID(instanceID),
     m_name(cdm_name){     
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
    if ( ! XmlUtil::writeIt(*doc,fileName) ) return 4;
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
  fsw_cdb_inst<SCHEMA>::fsw_cdb_inst(const char* cdm_name,
				     unsigned int key, short schemaID, 
				     short versionID, short instanceID, SCHEMA* value) 
    :fsw_cdb(cdm_name,key,schemaID,versionID,instanceID),
     m_schema(value){
  }
    
  template <typename SCHEMA>
  fsw_cdb_inst<SCHEMA>::~fsw_cdb_inst() {
    delete m_schema;
  }
    
  template <typename SCHEMA>
  fsw_datum* fsw_cdb_inst<SCHEMA>::get_io_handler() const {
    if ( m_schema == 0 ) return 0;
    std::string ioName;
    switch ( get_schemaID() ) {
    case fsw_enums::GFC_SCHEMA_ID:
      ioName = "GFC_DB";
      break;
    case fsw_enums::DFC_SCHEMA_ID:
      ioName = "DFC_DB";
      break;
    case fsw_enums::MFC_SCHEMA_ID:
      ioName = "MFC_DB"; 
      break;
    case fsw_enums::HFC_SCHEMA_ID:
      ioName = "HFC_DB"; 
      break;
    default:
      break;
    }
    fsw_datum* retVal = new fsw_datum_inst<SCHEMA>( ioName.c_str(), *m_schema );
    return retVal;
  }

}

