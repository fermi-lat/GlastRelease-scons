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
#include "./fsw_cdb_print.h"

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
  const SCHEMA* fsw_cdb_inst<SCHEMA>::load( int option ) {    
    const void* schema_ptr = fsw_get_data( get_name().c_str(), option,
					   &m_key, &m_schemaID, &m_versionID, &m_instanceID );
    m_schema = static_cast<const SCHEMA*>(schema_ptr);
    return m_schema;
  }
  
  template <typename SCHEMA>
  int fsw_cdb_inst<SCHEMA>::printToXmlFile( const char* fileName ) const {
    // print to xml file
    // return 0 for success, error flag otherwise
    //std::ofstream ofile(fileName);    
    int ok = genericPrint(std::cout);
    ok = specializedPrint(std::cout);
    //ofile.close();
    return ok;
  }

  template <typename SCHEMA>
  int fsw_cdb_inst<SCHEMA>::genericPrint( std::ostream& os ) const {
    os << std::endl;
    os << std::endl;
    os << "File: " << get_name() << std::endl;
    os << "CDM identification" << std::endl;
    os << " Key                 :  " << std::hex << m_key << std::endl;    
    os << " Schema              :  " << std::hex << m_schemaID << std::endl;
    os << " Version             : " << std::dec << m_versionID << std::endl;
    os << " Instance            : " << std::dec << m_instanceID << std::endl;
    return 0;
  }
 
  template <>
  int fsw_cdb_inst<GFC_DB_schema>::specializedPrint( std::ostream& os ) const {
    return fsw_cdb_print::print_gfc( get_schema(), os );
  }
  
  template <>
  int fsw_cdb_inst<DFC_DB_schema>::specializedPrint( std::ostream& os ) const {
    return fsw_cdb_print::print_dfc( get_schema(), os );
  } 
  
  template <>
  int fsw_cdb_inst<MFC_DB_schema>::specializedPrint( std::ostream& os ) const {
    return fsw_cdb_print::print_mfc( get_schema(), os );
  }
  
  template <>
  int fsw_cdb_inst<HFC_DB_schema>::specializedPrint( std::ostream& os ) const {
    return fsw_cdb_print::print_hfc( get_schema(), os );
  }
  
}
