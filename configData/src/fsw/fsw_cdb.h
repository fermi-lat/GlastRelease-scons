// -*- Mode: c++ -*-
#ifndef fsw_cdb_HH
#define fsw_cdb_HH
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Base class to wrap FSW CDM data and interface to the IO handlers
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



// c++/stl
#include <string>

// local includes

// forward declares
namespace configData {
  class fsw_datum;
}

namespace configData {

  struct fmx_id {
  public:    
    fmx_id(unsigned int key, short schemaID, 
	   short versionID, short instanceID)
      :m_key(key),m_schemaID(schemaID),
       m_versionID(versionID),m_instanceID(instanceID){
    }
    ~fmx_id(){;}    
    /// get the The FMX key of this CDM
    inline unsigned int get_key() const { return m_key; }
    /// get the ID of the schema of the CDM (ie, what kind of schema)
    inline short get_schemaID() const { return m_schemaID; }
    /// get the ID of the version of schema used for this CDM
    inline short get_versionID() const { return m_versionID; }
    /// get the ID of the instance of the CDM (relative to other CDM of same type)
    inline short get_instanceID() const { return m_instanceID; }
  private:
    /// The FMX key of this CDM
    unsigned int   m_key;     
    /// ID of the schema of the CDM (ie, what kind of schema)
    short          m_schemaID;   
    /// ID of the version of schema used for this CDM
    short          m_versionID;  
    /// ID of the instance of the CDM (relative to other CDM of same type)
    short          m_instanceID; 
  };


  /**  
   *
   *   @class fsw_cdb
   *   @author E. Charles
   *
   *   Base class for wrapping FSW CDMs with IO functionality
   *
   *   Stores bookkeeping information for any loaded CDM
   *   Also provides 3 IO functions
   *
   *      printToStream   -> text dump to a std::ostream
   *      writeToXmlFile  -> dump to xml
   *      readFromXmlFile -> overwrite values from an input xml file
   *
   **/

  class fsw_cdb {
    
  public:
    
    /**
     * @brief C'tor for fsw_cdb
     *
     * @param cdm_name is the file name of the shared object file
     * @param key is the FMX key of the CDM 
     * @param schemaID identifies the type of CDM
     * @param versionID identifies the version of that particular type of CDM
     * @param instanceID identifies the particular instance of CDM      
     **/   
    fsw_cdb(const char* cdm_name, 
	    unsigned int key, short schemaID, 
	    short versionID, short instanceID);

    /// D'tor for fsw_cdb
    virtual ~fsw_cdb();
    
    
    // access ----------------------------------
    
    /// get the name of the CDM
    inline const std::string& get_name() const { return m_name; }
    /// get the FMX information
    inline const fmx_id& get_fmx_id() const { return m_fmxId; }

    /**
     * @brief print data to a stream
     *
     * @param os is the stream
     * @return 0 for success, failure code otherwise
     **/   
    int printToStream( std::ostream& os ) const;

    /**
     * @brief write to an XML file.  overwrites existing files
     *
     * @param fileName is the name of the file
     * @return 0 for success, failure code otherwise
     **/   
    int writeToXmlFile( const char* fileName ) const;
   
    /**
     * @brief overwrite values in memory with values from and xml file
     *
     * @param fileName is the name of the file
     * @return 0 for success, failure code otherwise
     **/   
    int readFromXmlFile( const char* fileName );
       
  protected:

    /// Make the object that handles the IO for this schema
    virtual fsw_datum* get_io_handler() const = 0;

  private: 
    
    /// The FMX bookkeeping information
    fmx_id         m_fmxId;

    /// The name of the CDM we are printing
    std::string    m_name;        
  }; 

}

#endif
