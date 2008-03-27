// -*- Mode: c++ -*-
#ifndef fsw_datum_inst_HH
#define fsw_datum_inst_HH
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Templated class for converting FSW headers to XML or text
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
#include "./fsw_datum.h"

// c++/stl
#include <string>
#include <list>
#include <map>
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


  /** 
   *  @class fsw_datum_inst
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
   *  For simple types the printout looks like
   *    var_name              : value
   *
   *  For structures the printout looks like
   *    struct_name
   *      var_name_1          : value1
   *      var_name_2          : value2
   *
   *  For simple types the XML looks like
   *    <var_name>value</var_name>
   *
   *  For structures the printout looks like
   *    <struct_name>
   *      <var_name_1>value1</var_name_1>
   *      <var_name_2>value2</var_name_2>
   *    </struct_name>
   *
   **/

  template <typename T>
  class fsw_datum_inst : public fsw_datum {
    
  public:
    
    /**
     * @brief C'tor for fsw_datum_inst, gives a name to some data
     *
     * This name will be used in both the XML and text printouts
     *
     * This class can handle data with sub-structure.
     *
     *   addChild() will add a sub-element
     *   getChildByName() will find a sub-element
     *
     * @param name is the name associated with the datum
     * @param val is a reference to the data which will be managed by this class
     **/
    fsw_datum_inst(const char* name, T& val);
    
    /// D'tor
    virtual ~fsw_datum_inst();
    
    
    // access ----------------------------------
    
    /// get the data
    inline const T& get_datum() const { return m_datum; }
    

    /**
     * @brief get a sub-element of the data structure
     *
     * @param name is the name of the sub-element
     * @return a pointer to the IO handler for the sub-element of data, null if no such child
     **/
    fsw_datum* getChildByName( const char* name );


    // modifiers -------------------------------

    /**
     * @brief add sub-element of the data structure
     *
     * @param datum is the IO handler for the sub-element
     **/
    void addChild( fsw_datum& datum );


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

    
  protected:    

    /// Print the children
    void printChildren( int depth, std::ostream& os = std::cout) const;
    
    /// Write the children to XML
    bool writeChildren( DOMElement& thisNode ) const;

    /// Read the children from XML
    bool readChildren( DOMElement& thisNode );

    /**
     * @brief print to a stream, varies depending on template specialization
     *
     * @param depth is the indent depth
     * @param os is the stream
     **/
    void specializedPrint( int depth, std::ostream& os = std::cout ) const;

    /**
     * @brief create and fill a XML element with this datum
     * varies depending on template specialization
     *
     * @param parent is the parent xml element
     * @return the newly created XML element
     **/
    DOMElement* specializedWriteToXml( DOMElement& parent ) const;


    /**
     * @brief set CDM values from an XML element
     * varies depending on template specialization
     *
     * @param node is the xml element
     * @return true for success, false otherwise
     **/
    bool specializedReadFromXml( DOMElement& thisNode );


  private:

    /// Reference to the structure handled by this object
    T&                               m_datum;    

    /// All the children of this datum
    std::list<fsw_datum*>            m_children;
    
    /// All the children, sorted by name
    std::map<std::string,fsw_datum*> m_childMap;


  };

  // Simple name for all the basic types we use

  // pointer
  typedef fsw_datum_inst<void*>                 fsw_datum_ptr;
  // 32 bit unsigned int
  typedef fsw_datum_inst<unsigned int>          fsw_datum_uint;

}


#ifndef fsw_datum_inst_CXX
#include "./fsw_datum_inst.cxx"
#endif

#endif
