// -*- Mode: c++ -*-
#ifndef fsw_cdb_HH
#define fsw_cdb_HH
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
/** @class fsw_cdb
    @author E. Charles

    .Base class for converting FSW headers to XML
*/


// c++/stl
#include <string>

// local includes
#include "./fsw_enums.h"
#include "./fsw_iface.h"


// forward declares
namespace configData {
  class fsw_datum;
}

namespace configData {

  class fsw_cdb {
    
  public:
    
    static fsw_cdb* read_cdb( fsw_enums::CDM_TYPE cdmType, const char* name, int option = 0 );
    
    static fsw_enums::CDM_TYPE cdm_type( const std::string& cdb_type_name );

  public:
    
    // C'tor
    fsw_cdb(const char* cdm_name);
    
    // D'tor
    virtual ~fsw_cdb();
    
    
    // access ----------------------------------
    
    // get the name of the CDM
    inline const std::string& get_name() const { return m_name; }
    
    // print to a stream
    // return 0 for success, error flag otherwise
    int printToStream( std::ostream& os ) const;

    // write to an XML file
    // return 0 for success, error flag otherwise
    int writeToXmlFile( const char* fileName ) const;
   
    // write from an XML file
    // return 0 for success, error flag otherwise    
    int readFromXmlFile( const char* fileName );
       
  protected:

    virtual fsw_datum* get_io_handler() const = 0;

    unsigned int m_key;
    short        m_schemaID;
    short        m_versionID;
    short        m_instanceID;

  private: 
    
    std::string    m_name;     // The name of the CDM we are printing
  }; 



  template <typename SCHEMA>
  class fsw_cdb_inst : public fsw_cdb {
    
  public:
    
    // C'tor
    fsw_cdb_inst(const char* cdm_name, int option = 0);
    
    // D'tor
    virtual ~fsw_cdb_inst();
    
    
    // access ----------------------------------
    
    // get the data
    inline const SCHEMA* get_schema() const { return m_schema; }
    // get the data
    inline SCHEMA* get_schema() { return m_schema; }
    
    
    // class methods
    
    // Read the data and cast it to the correct SCHEMA
    // null pointer means failure
    SCHEMA* load( int option = 0 );
       
  protected:

    virtual fsw_datum* get_io_handler() const;

  private:
    
    mutable SCHEMA*  m_schema;    
  };

}

#endif
