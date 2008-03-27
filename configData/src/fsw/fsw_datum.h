// -*- Mode: c++ -*-
#ifndef fsw_datum_HH
#define fsw_datum_HH
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
/** @class fsw_datum
    @author E. Charles

    .Base class for converting FSW headers to XML
*/



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


  // Base class for managing IO on FSW CDM

  class fsw_datum {
   
  public:
    
    // C'tor
    fsw_datum(const char* name);
    
    // D'tor
    virtual ~fsw_datum();
    
    
    // access ----------------------------------
    
    // get the name of the CDM
    inline const std::string& get_name() const { return m_name; }
  

    // virtual methods
    
    // print to a stream
    virtual void print( int depth, std::ostream& os = std::cout ) const = 0;

    // Write to XML
    virtual DOMElement* writeToXml( DOMElement& parent ) const = 0;

    // Read from XML
    virtual bool readFromXml( DOMElement& node ) = 0;
    
  protected:

  private: 

    // The name of this datum
    std::string m_name;     
  }; 




  class fsw_datum_mask : public fsw_datum {

  public:
    
    // C'tor
    fsw_datum_mask(const char* name, unsigned int& val);
    
    // D'tor
    virtual ~fsw_datum_mask();

    // access ----------------------------------
    
    // get the data
    inline unsigned int get_mask() const { return m_mask; }

    // print to a stream
    virtual void print( int depth, std::ostream& os = std::cout ) const;

    // Write to XML
    virtual DOMElement* writeToXml( DOMElement& parent ) const;

    // Read from XML
    virtual bool readFromXml( DOMElement& node );

  private:

    unsigned int& m_mask;
  };


  template <typename T>
  class fsw_datum_inst : public fsw_datum {
    
  public:
    
    // C'tor
    fsw_datum_inst(const char* name, T& val);
    
    // D'tor
    virtual ~fsw_datum_inst();
    
    
    // access ----------------------------------
    
    // get the data
    inline const T& get_datum() const { return m_datum; }
    

    // get a child
    fsw_datum* getChildByName( const char* name );


    // modifiers -------------------------------

    // add a child node
    void addChild( fsw_datum& datum );


    // print to a stream
    virtual void print( int depth, std::ostream& os = std::cout ) const;

    // Write to XML
    virtual DOMElement* writeToXml( DOMElement& parent ) const;

    // Read from XML
    virtual bool readFromXml( DOMElement& node );
    
  protected:
    

    // Print the children
    void printChildren( int depth, std::ostream& os = std::cout) const;
    
    // Write the children
    bool writeChildren( DOMElement& thisNode ) const;

    // Read the children
    bool readChildren( DOMElement& thisNode );

    // print to a stream
    void specializedPrint( int depth, std::ostream& os = std::cout ) const;

    // Write to XML
    DOMElement* specializedWriteToXml( DOMElement& parent ) const;

    // Read from XML
    bool specializedReadFromXml( DOMElement& thisNode );


  private:

    // This datum
    T&                               m_datum;    

    // All the children of this datum
    std::list<fsw_datum*>            m_children;
    
    // All the children, sorted by name
    std::map<std::string,fsw_datum*> m_childMap;


  };

  // All the basic types we use
  typedef fsw_datum_inst<void*>                 fsw_datum_ptr;
  typedef fsw_datum_inst<unsigned int>          fsw_datum_uint;
  typedef fsw_datum_inst<unsigned short>        fsw_datum_ushort;  

}

#endif
