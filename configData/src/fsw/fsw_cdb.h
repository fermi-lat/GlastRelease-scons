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
    
    // This is the important function
    // 
    // Load a CDM, cast to the correct type and return us 
    // the base class
    //
    // cdm_name is the file name of the shared object file
    // option is not used 
    static fsw_cdb* read_cdb( const char* cdm_name, int option = 0 );
    
  public:
    
    // C'tor
    fsw_cdb(const char* cdm_name, 
	    unsigned int key, short schemaID, 
	    short versionID, short instanceID);
    
    // D'tor
    virtual ~fsw_cdb();
    
    
    // access ----------------------------------
    
    // get the name of the CDM
    inline const std::string& get_name() const { return m_name; }
    // get the The FMX key of this CDM
    inline unsigned int get_key() const { return m_key; }
    // get the ID of the schema of the CDM (ie, what kind of schema)
    inline short get_schemaID() const { return m_schemaID; }
    // get the ID of the version of schema used for this CDM
    inline short get_versionID() const { return m_versionID; }
    // get the ID of the instance of the CDM (relative to other CDM of same type)
    inline short get_instanceID() const { return m_instanceID; }


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

    // Make the object that handles the IO for this schema
    virtual fsw_datum* get_io_handler() const = 0;

  private: 
    
    // The FMX key of this CDM
    unsigned int   m_key;     
    // ID of the schema of the CDM (ie, what kind of schema)
    short          m_schemaID;   
    // ID of the version of schema used for this CDM
    short          m_versionID;  
    // ID of the instance of the CDM (relative to other CDM of same type)
    short          m_instanceID; 
    // The name of the CDM we are printing
    std::string    m_name;        
  }; 



  template <typename SCHEMA>
  class fsw_cdb_inst : public fsw_cdb {
    
  public:
    
    // C'tor
    fsw_cdb_inst(const char* cdm_name,  
		 unsigned int key, short schemaID, 
		 short versionID, short instanceID,
		 SCHEMA* value);
    
    // D'tor
    virtual ~fsw_cdb_inst();
    
    
    // access ----------------------------------
    
    // get the data
    inline const SCHEMA* get_schema() const { return m_schema; }
    // get the data
    inline SCHEMA* get_schema() { return m_schema; }

              
  protected:

    // Make the object that handles the IO for this schema
    virtual fsw_datum* get_io_handler() const;

  private:
    
    // The data
    mutable SCHEMA*  m_schema;    
  };

}

#endif
