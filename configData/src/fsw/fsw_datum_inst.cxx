#ifndef fsw_datum_inst_CXX
#define fsw_datum_inst_CXX
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
#include "./fsw_datum_inst.h"

// c++/stl headers
#include <iostream>
#include <iomanip>

// Local headers
#include "configData/base/ConfigXmlUtil.h"



namespace configData {

  // C'tor for fsw_datum_inst, gives a name to some data
  template <typename T>
  fsw_datum_inst<T>::fsw_datum_inst(const char* name, T& val)
    :fsw_datum(name),
     m_datum(val){
  }

  // D'tor must delete children
  template <typename T>
  fsw_datum_inst<T>::~fsw_datum_inst() {
    for ( std::list<fsw_datum*>::iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      delete child;
    }
  }

  // Get a child data element by name, return 0 for failure
  template <typename T>
  fsw_datum* fsw_datum_inst<T>::getChildByName( const char* name ) {
    std::map<std::string,fsw_datum*>::iterator itr = m_childMap.find(name);
    return itr != m_childMap.end() ? itr->second : 0;
  }

  // Add a child data element
  template <typename T>
  void fsw_datum_inst<T>::addChild( fsw_datum& datum ) {
    // add to the list of childern
    m_children.push_back( &datum );
    // add to the map, by name
    m_childMap[datum.get_name()] = &datum;
  }

  // print to a stream
  template <typename T>
  void fsw_datum_inst<T>::print( int depth, std::ostream& os ) const {
    // print self
    specializedPrint(depth, os);
    // print children
    printChildren(depth, os);
  }

  // create and fill a XML element with this datum
  template <typename T>
  DOMElement* fsw_datum_inst<T>::writeToXml( DOMElement& parent ) const {
    // create element for self
    DOMElement* elem = specializedWriteToXml(parent);
    if ( elem == 0 ) return 0;
    // create sub-elements for children
    if ( ! writeChildren(*elem) ) return 0;
    return elem;
  }
  
  // create and fill a XML element with this datum
  template <typename T>
  bool fsw_datum_inst<T>::readFromXml( DOMElement& node ) {
    // fill values for self
    if ( ! specializedReadFromXml(node) ) return false;
    // fill values for children
    if ( ! readChildren(node) ) return false;
    return true;
  }

  // print children to a std::ostream
  template <typename T>
  void fsw_datum_inst<T>::printChildren( int depth, std::ostream& os ) const {
    // up the indent level
    int newDepth = depth + 2;
    // loop on children and print them
    for ( std::list<fsw_datum*>::const_iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      child->print(newDepth,os);
    }
  }

  // create XML sub-elements for children
  template <typename T>
  bool fsw_datum_inst<T>::writeChildren( DOMElement& thisNode ) const {
    // loop on children
    for ( std::list<fsw_datum*>::const_iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      // make sub-elements
      DOMElement* ok = child->writeToXml(thisNode);
      if ( ok == 0 ) return false;
    }
    return true;
  }

  // fill values for children
  template <typename T>
  bool fsw_datum_inst<T>::readChildren( DOMElement& thisNode ) {
    // loop on children
    for ( std::map<std::string,fsw_datum*>::iterator itr = m_childMap.begin();
	  itr != m_childMap.end(); itr++ ) {
      // get the right XML node
      DOMElement* childNode = XmlUtil::findChildByName( thisNode, itr->first.c_str() );
      if ( childNode == 0 ) return false;
      fsw_datum* child = itr->second;
      // fill values for child element from XML
      bool ok = child->readFromXml(*childNode);
      if ( !ok ) return false;
    }
    return true;
  }
  
  // Default version of print to a stream
  template <typename T>
  void fsw_datum_inst<T>::specializedPrint( int depth, std::ostream& os ) const {
    // just print the element name
    os << std::setw(depth) << ' ' << get_name() << std::endl;
  }
  
  // Default version of write to XML
  template <typename T>
  DOMElement* fsw_datum_inst<T>::specializedWriteToXml( DOMElement& parent ) const {
    // just make a node with the element name
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );    
    return thisNode;
  }

  //  // Default version of read from XML
  template <typename T>
  bool fsw_datum_inst<T>::specializedReadFromXml( DOMElement& thisNode ) {
    // Already have the node name.  done
    return true;
  }

  // declare partial specializations for simple types

  // Print a pointer to a stream
  template <>
  void fsw_datum_inst<void*>::specializedPrint( int depth, std::ostream& os) const;

  // Print an unsigned int to a stream
  template <>
  void fsw_datum_inst<unsigned int>::specializedPrint( int depth, std::ostream& os) const;

  // Write a pointer xml
  template <>
  DOMElement* fsw_datum_inst<void*>::specializedWriteToXml( DOMElement& parent ) const;

  // Write an unsigned int to xml
  template <>
  DOMElement* fsw_datum_inst<unsigned int>::specializedWriteToXml( DOMElement& parent ) const;
  
  // Read a pointer from XML
  template <>
  bool fsw_datum_inst<void*>::specializedReadFromXml( DOMElement& thisNode );
  
  // Read an unsigned int from XML
  template <>
  bool fsw_datum_inst<unsigned int>::specializedReadFromXml( DOMElement& thisNode );

}


// forward declare the various schema
struct _EFC_DB_prescalers;
struct _EFC_DB_sampler;

struct _DFC_DB_gemV0;
struct _DFC_DB_stats;
struct _DFC_DB_prms;
struct _DFC_DB_schema;

struct _MFC_DB_calV0;
struct _MFC_DB_tkrV0;
struct _MFC_DB_prms;
struct _MFC_DB_schema;

struct _HFC_DB_acdV0;
struct _HFC_DB_calV0;
struct _HFC_DB_tkrV0;
struct _HFC_DB_prms;
struct _HFC_DB_schema;

struct _GFC_DB_acdV0;
struct _GFC_DB_atfV0;
struct _GFC_DB_zbottomV0;
struct _GFC_DB_calV1;
struct _GFC_DB_tkrV0;
struct _GFC_DB_prms;
struct _GFC_DB_schema;


// Forward declare the partial specializations for the schema
namespace configData {

  // C'tors that add children
  template <>
  fsw_datum_inst<_EFC_DB_prescalers>::fsw_datum_inst(const char* name, _EFC_DB_prescalers& val);

  template <>
  fsw_datum_inst<_EFC_DB_sampler>::fsw_datum_inst(const char* name, _EFC_DB_sampler& val);

  template <>
  fsw_datum_inst<_DFC_DB_gemV0>::fsw_datum_inst(const char* name, _DFC_DB_gemV0& val);
 
  template <>
  fsw_datum_inst<_DFC_DB_stats>::fsw_datum_inst(const char* name, _DFC_DB_stats& val);

  template <>
  fsw_datum_inst<_DFC_DB_prms>::fsw_datum_inst(const char* name, _DFC_DB_prms& val);

  template <>
  fsw_datum_inst<_MFC_DB_calV0>::fsw_datum_inst(const char* name, _MFC_DB_calV0& val);

  template <>
  fsw_datum_inst<_MFC_DB_tkrV0>::fsw_datum_inst(const char* name, _MFC_DB_tkrV0& val);

  template <>
  fsw_datum_inst<_MFC_DB_prms>::fsw_datum_inst(const char* name, _MFC_DB_prms& val);

  template <>
  fsw_datum_inst<_HFC_DB_acdV0>::fsw_datum_inst(const char* name, _HFC_DB_acdV0& val);

  template <>
  fsw_datum_inst<_HFC_DB_calV0>::fsw_datum_inst(const char* name, _HFC_DB_calV0& val);

  template <>
  fsw_datum_inst<_HFC_DB_tkrV0>::fsw_datum_inst(const char* name, _HFC_DB_tkrV0& val);

  template <>
  fsw_datum_inst<_HFC_DB_prms>::fsw_datum_inst(const char* name, _HFC_DB_prms& val);

  template <>
  fsw_datum_inst<_GFC_DB_acdV0>::fsw_datum_inst(const char* name, _GFC_DB_acdV0& val);

  template <>
  fsw_datum_inst<_GFC_DB_atfV0>::fsw_datum_inst(const char* name, _GFC_DB_atfV0& val);

  template <>
  fsw_datum_inst<_GFC_DB_zbottomV0>::fsw_datum_inst(const char* name, _GFC_DB_zbottomV0& val);

  template <>
  fsw_datum_inst<_GFC_DB_calV1>::fsw_datum_inst(const char* name, _GFC_DB_calV1& val);

  template <>
  fsw_datum_inst<_GFC_DB_tkrV0>::fsw_datum_inst(const char* name, _GFC_DB_tkrV0& val);

  template <>
  fsw_datum_inst<_GFC_DB_prms>::fsw_datum_inst(const char* name, _GFC_DB_prms& val);

  template <>
  fsw_datum_inst<_GFC_DB_schema>::fsw_datum_inst(const char* name, _GFC_DB_schema& val);

  template <>
  fsw_datum_inst<_HFC_DB_schema>::fsw_datum_inst(const char* name, _HFC_DB_schema& val);

  template <>
  fsw_datum_inst<_MFC_DB_schema>::fsw_datum_inst(const char* name, _MFC_DB_schema& val);

  template <>
  fsw_datum_inst<_DFC_DB_schema>::fsw_datum_inst(const char* name, _DFC_DB_schema& val);


  // Specialized printouts
  template <>
  void fsw_datum_inst<_EFC_DB_prescalers>::specializedPrint( int depth, std::ostream& os) const;
   
  template <>
  void fsw_datum_inst<_DFC_DB_gemV0>::specializedPrint( int depth, std::ostream& os) const;

  template <>
  void fsw_datum_inst<_DFC_DB_stats>::specializedPrint( int depth, std::ostream& os) const;


  // Specialized XML formatting
  template <>
  DOMElement* fsw_datum_inst<_EFC_DB_prescalers>::specializedWriteToXml( DOMElement& parent ) const;

  template <>
  DOMElement* fsw_datum_inst<_DFC_DB_gemV0>::specializedWriteToXml( DOMElement& parent ) const;

  template <>
  DOMElement* fsw_datum_inst<_DFC_DB_stats>::specializedWriteToXml( DOMElement& parent ) const;


  // Specialized XML reading
  template <>
  bool fsw_datum_inst<_EFC_DB_prescalers>::specializedReadFromXml( DOMElement& thisNode );

  template <>
  bool fsw_datum_inst<_DFC_DB_gemV0>::specializedReadFromXml( DOMElement& thisNode );

  template <>
  bool fsw_datum_inst<_DFC_DB_stats>::specializedReadFromXml( DOMElement& thisNode );

}

#endif
