// -*- Mode: c++ -*-
#ifndef fsw_cdb_inst_HH
#define fsw_cdb_inst_HH
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Templated class to wrap FSW CDM data and interface to the IO handlers
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

// base class header
#include "./fsw_cdb.h"

// c++/stl
#include <string>

// local includes
#include "./fsw_iface.h"


// forward declares
namespace configData {
  class fsw_datum;
}

namespace configData {

  /**  
   *
   *   @class fsw_cdb_inst
   *   @author E. Charles
   *
   *   Template for handling specific CDMs
   *   This is here to cast the CDM struct to the correct type
   *   And to correctly implement the get_io_handler() function
   *
   **/

  template <typename SCHEMA>
  class fsw_cdb_inst : public fsw_cdb {
    
  public:
    
    /**
     * @brief C'tor for fsw_cdb_inst
     *
     * @param cdm_name is the file name of the shared object file
     * @param key is the FMX key of the CDM 
     * @param schemaID identifies the type of CDM
     * @param versionID identifies the version of that particular type of CDM
     * @param instanceID identifies the particular instance of CDM           
     * @param value points to the memory where the CDM has been loaded           
     **/   
    fsw_cdb_inst(const char* cdm_name,  
		 unsigned int key, short schemaID, 
		 short versionID, short instanceID,
		 SCHEMA* value);
    
    /// D'tor for fsw_cdb_inst
    virtual ~fsw_cdb_inst();
    
    
    // access ----------------------------------
    
    /// get the CDM data
    inline const SCHEMA* get_schema() const { return m_schema; }
    /// get the CDM data
    inline SCHEMA* get_schema() { return m_schema; }

              
  protected:

    /// Make the object that handles the IO for this schema
    virtual fsw_datum* get_io_handler() const;

  private:
    
    /// Pointer to the CDM data, correctly cast
    mutable SCHEMA*  m_schema;    
  };

}

#ifndef fsw_cdb_inst_CXX
#include "./fsw_cdb_inst.cxx"
#endif

#endif
