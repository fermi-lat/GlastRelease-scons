// -*- Mode: c++ -*-
#ifndef fsw_datum_HH
#define fsw_datum_HH
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Base class for converting FSW headers to XML or text
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
#include <iostream>

// stuff for xerces
#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class  DOMElement;
XERCES_CPP_NAMESPACE_END

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

// local includes

// forward declares
namespace configData {
  class fsw_cdb;
  struct fmx_id;
}

namespace configData {


  /** 
   *  @class fsw_datum
   *  @author E. Charles
   *
   *  Wrap data for IO functionality 
   *  The idea is to wrap any structure used in CDM
   *
   *  This class define 3 virtual funtions
   *
   *   print        -> to a std::ostream
   *   writeToXml   -> make and return a DOMElement
   *   readFromXml  -> set values from a DOMElement
   *
   *  At the lowest level this just wraps simple types like unsigned int
   *  At higher levels this wraps sub-structure of the CDM
   *  At the top level this wraps the CDM itself
   *
   **/

  class fsw_datum {

  public:    
  
    /**
     * @brief Load a CDM, cast to the correct type and return us 
     * the base class
     * 
     * @param cdm_name is the file name of the shared object file
     * @param option is not used 
     * @return an fsw_cdb that wraps the CDM, null for failure
     **/
    static fsw_cdb* read_cdb( const char* cdm_name, int option = 0 );
   
    /**
     * @brief write bookkeeping information to an XML file
     * 
     * @param parent is the parent xml element
     * @param cdm_name is the file name of the shared object file
     * @param fmxId is the FMX bookkeeping information 
     * @return the newly created XML element
     **/
    static DOMElement* writeXmlHeader( DOMElement& parent, const char* cdm_name, const fmx_id& fmxId );

  public:

    /**
     * @brief C'tor for fsw_datum gives a name to any item of data
     *
     * This name will be used in both the XML and text printouts
     *
     * @param name is the name associated with the datum
     **/
    fsw_datum(const char* name);
    
    /// D'tor for fsw_datum
    virtual ~fsw_datum();
    
    
    // access ----------------------------------
    
    /// get the name this datum
    inline const std::string& get_name() const { return m_name; }
  

    // virtual methods
    
    /**
     * @brief print to a stream
     *
     * @param depth is the indent depth
     * @param os is the stream
     **/
    virtual void print( int depth, std::ostream& os = std::cout ) const = 0;

    /**
     * @brief create and fill a XML element with this datum
     *
     * @param parent is the parent xml element
     * @return the newly created XML element
     **/
    virtual DOMElement* writeToXml( DOMElement& parent ) const = 0;

    /**
     * @brief set CDM values from an XML element
     *
     * @param node is the xml element
     * @return true for success, false otherwise
     **/
    virtual bool readFromXml( DOMElement& node ) = 0;
    
  protected:

  private: 

    /// The name of this datum
    std::string m_name;     
  }; 


  /** 
   *  @class fsw_datum_mask
   *  @author E. Charles
   *
   *  Wrap 32 bit unsigned int for IO
   *
   *  This class implements the io functions
   *
   *   print        -> to a std::ostream as a series of '.' (0) or 'X' (1)
   *      name                    :.....XXXX...XX.....XX...XX.X.X.X
   *   writeToXml   -> make and return a DOMElement
   *   readFromXml  -> set values from a DOMElement
   *      <name>0x5ff23091</name>
   *
   **/

  class fsw_datum_mask : public fsw_datum {

  public:
    
    /**
     * @brief C'tor for fsw_datum_mask, gives a name to a 32 bit mask
     *
     * This name will be used in both the XML and text printouts
     *
     * @param name is the name associated with the datum
     * @param val is a reference to the mask, will be managed by this class
     **/
    fsw_datum_mask(const char* name, unsigned int& val);
    
    /// D'tor
    virtual ~fsw_datum_mask();


    // access ----------------------------------
    
    /// get the data
    inline unsigned int get_mask() const { return m_mask; }


    // virtual methods -------------------------

    /**
     * @brief print to a stream
     *
     * @param depth is the indent depth
     * @param os is the stream
     **/
    virtual void print( int depth, std::ostream& os = std::cout ) const;

    
    /**
     * @brief create and fill a XML element with this datum
     *
     * @param parent is the parent xml element
     * @return the newly created XML element
     **/
    virtual DOMElement* writeToXml( DOMElement& parent ) const;

    /**
     * @brief set CDM values from an XML element
     *
     * @param node is the xml element
     * @return true for success, false otherwise
     **/
    virtual bool readFromXml( DOMElement& node );

  private:

    /// reference to 32 bit mask handled by this object
    unsigned int& m_mask;
  };

}

#endif
