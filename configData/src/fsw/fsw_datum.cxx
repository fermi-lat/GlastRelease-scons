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
#include "./fsw_datum.h"

// c++/stl headers
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>


// Local headers
#include "configData/base/ConfigXmlUtil.h"



namespace configData {

  fsw_datum::fsw_datum(const char* name) 
    :m_name(name){
  }
  
  fsw_datum::~fsw_datum() {
  }
  
  

  fsw_datum_mask::fsw_datum_mask(const char* name, unsigned int& val)
    :fsw_datum(name),
     m_mask(val){
  }

  fsw_datum_mask::~fsw_datum_mask() {
  }


  void fsw_datum_mask::print( int depth, std::ostream& os ) const {
    os << std::setw(depth) << ' '  
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ": ";
    unsigned int test(0x1);
    for ( int i(0); i < 32; i++ ) {
      os << ( m_mask & test ? '.' : 'X' );
      test *= 2;
    }
    os << std::endl;
  }

  DOMElement* fsw_datum_mask::writeToXml( DOMElement& parent ) const {
    char tmp[12];
    sprintf(tmp,"0x%8x",m_mask);
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), tmp);
  }
  
  bool fsw_datum_mask::readFromXml( DOMElement& node ) {
    return true;
  }
  

  template <typename T>
  fsw_datum_inst<T>::fsw_datum_inst(const char* name, T& val)
    :fsw_datum(name),
     m_datum(val){
  }

  template <typename T>
  fsw_datum_inst<T>::~fsw_datum_inst() {
    for ( std::list<fsw_datum*>::iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      delete child;
    }
  }

  template <typename T>
  fsw_datum* fsw_datum_inst<T>::getChildByName( const char* name ) {
    std::map<std::string,fsw_datum*>::iterator itr = m_childMap.find(name);
    return itr != m_childMap.end() ? itr->second : 0;
  }

  template <typename T>
  void fsw_datum_inst<T>::addChild( fsw_datum& datum ) {
    m_children.push_back( &datum );
    m_childMap[datum.get_name()] = &datum;
  }

  template <typename T>
  void fsw_datum_inst<T>::print( int depth, std::ostream& os ) const {
    specializedPrint(depth, os);
    printChildren(depth, os);
  }

  template <typename T>
  DOMElement* fsw_datum_inst<T>::writeToXml( DOMElement& parent ) const {
    DOMElement* elem = specializedWriteToXml(parent);
    if ( elem == 0 ) return 0;
    if ( ! writeChildren(*elem) ) return 0;
    return elem;
  }
  
  template <typename T>
  bool fsw_datum_inst<T>::readFromXml( DOMElement& node ) {
    if ( ! specializedReadFromXml(node) ) return false;
    if ( ! readChildren(node) ) return false;
    return true;
  }

  template <typename T>
  void fsw_datum_inst<T>::printChildren( int depth, std::ostream& os ) const {
    int newDepth = depth + 2;
    for ( std::list<fsw_datum*>::const_iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      child->print(newDepth,os);
    }
  }

  template <typename T>
  bool fsw_datum_inst<T>::writeChildren( DOMElement& thisNode ) const {
    for ( std::list<fsw_datum*>::const_iterator itr = m_children.begin();
	  itr != m_children.end(); itr++ ) {
      fsw_datum* child = *itr;
      DOMElement* ok = child->writeToXml(thisNode);
      if ( ok == 0 ) return false;
    }
    return true;
  }

  template <typename T>
  bool fsw_datum_inst<T>::readChildren( DOMElement& thisNode ) {
    for ( std::map<std::string,fsw_datum*>::iterator itr = m_childMap.begin();
	  itr != m_childMap.end(); itr++ ) {
      DOMElement* childNode = XmlUtil::findChildByName( thisNode, itr->first.c_str() );
      if ( childNode == 0 ) return false;
      fsw_datum* child = itr->second;
      bool ok = child->readFromXml(*childNode);
      if ( !ok ) return false;
    }
    return true;
  }
  
  // print to a stream
  template <typename T>
  void fsw_datum_inst<T>::specializedPrint( int depth, std::ostream& os ) const {
    os << std::setw(depth) << ' ' << get_name() << std::endl;
  }
  
  // Write to XML
  template <typename T>
  DOMElement* fsw_datum_inst<T>::specializedWriteToXml( DOMElement& parent ) const {
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );    
    return thisNode;
  }

  // Read from XML
  template <typename T>
  bool fsw_datum_inst<T>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

  // Print
  template <>
  void fsw_datum_inst<void*>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' '  
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ':' << m_datum << std::endl;
  }

  template <>
  void fsw_datum_inst<unsigned int>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' ' 
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ':' << m_datum << std::endl;
  }

  template <>
  void fsw_datum_inst<unsigned short>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' ' 
       << std::setw(35-depth) << std::left << get_name() << std::right
       << ':' << m_datum << std::endl;
  }

  // Write to XML
  template <>
  DOMElement* fsw_datum_inst<void*>::specializedWriteToXml( DOMElement& parent ) const {
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), m_datum);
  }

  template <>
  DOMElement* fsw_datum_inst<unsigned int>::specializedWriteToXml( DOMElement& parent ) const {
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), m_datum );
  }

  template <>
  DOMElement* fsw_datum_inst<unsigned short>::specializedWriteToXml( DOMElement& parent ) const {
    return XmlUtil::makeChildNodeWithContent( parent, get_name().c_str(), m_datum );
  }

  // Read from XML
  template <>
  bool fsw_datum_inst<void*>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

  template <>
  bool fsw_datum_inst<unsigned int>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

  template <>
  bool fsw_datum_inst<unsigned short>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

}

#include "./fsw_schema.icxx"
