#ifndef fsw_cdb_inst_CXX
#define fsw_cdb_inst_CXX
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Templated class for converting FSW headers to XML
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
#include "./fsw_cdb_inst.h"

// c++/stl headers
#include <iostream>
#include <iomanip>
#include <fstream>

// other headers
#include "configData/base/ConfigXmlUtil.h"
#include "./fsw_datum_inst.h"
#include "./fsw_enums.h"


namespace configData {

  // C'tor for fsw_cdb_inst
  template <typename SCHEMA>
  fsw_cdb_inst<SCHEMA>::fsw_cdb_inst(const char* cdm_name,
				     unsigned int key, short schemaID, 
				     short versionID, short instanceID, SCHEMA* value) 
    :fsw_cdb(cdm_name,key,schemaID,versionID,instanceID),
     m_schema(value){
  }
    
  // D'tor for fsw_cdb_inst
  template <typename SCHEMA>
  fsw_cdb_inst<SCHEMA>::~fsw_cdb_inst() {
    delete m_schema;
  }
   

  /// Make the object that handles the IO for this schema
  template <typename SCHEMA>
  fsw_datum* fsw_cdb_inst<SCHEMA>::get_io_handler() const {    
    // make sure we have some data to handle
    if ( m_schema == 0 ) return 0;
    // Get the name of the schema
    std::string ioName;
    switch ( get_fmx_id().get_schemaID() ) {
    case fsw_enums::GFC_SCHEMA_ID:
      ioName = "GFC_DB_schema";
      break;
    case fsw_enums::DFC_SCHEMA_ID:
      ioName = "DFC_DB_schema";
      break;
    case fsw_enums::MFC_SCHEMA_ID:
      ioName = "MFC_DB_schema"; 
      break;
    case fsw_enums::HFC_SCHEMA_ID:
      ioName = "HFC_DB_schema"; 
      break;
    default:
      break;
    }
    // make a handler for the type of schema we have, with the right name
    fsw_datum* retVal = new fsw_datum_inst<SCHEMA>( ioName.c_str(), *m_schema );
    return retVal;
  }

}

#endif
