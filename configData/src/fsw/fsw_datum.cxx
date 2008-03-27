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

// FSW headers
#include "GFC_DB/GFC_DB_schema.h"
#include "XFC_DB/HFC_DB_schema.h"
#include "XFC_DB/MFC_DB_schema.h"
#include "XFC_DB/DFC_DB_schema.h"


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


  // Specializations  
  
  template <>
  fsw_datum_inst<EFC_DB_prescalers>::fsw_datum_inst(const char* name, EFC_DB_prescalers& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Input_Prescale",val.in) ) );
    addChild( *(new fsw_datum_uint("Output_Prescale",val.out) ) );
    //addChild( *(new fsw_datum_uint_vect_32("Prescales",val.pv) ) );
  }

  template <>
  fsw_datum_inst<EFC_DB_sampler>::fsw_datum_inst(const char* name, EFC_DB_sampler& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_mask("Enabled",val.enabled) ) );
    addChild( *(new fsw_datum_inst<EFC_DB_prescalers>("Prescalers",val.prescalers) ) );
  }
    
  template <>
  fsw_datum_inst<DFC_DB_gemV0>::fsw_datum_inst(const char* name, DFC_DB_gemV0& val)
    :fsw_datum(name),
     m_datum(val){
    //addChild( *(new fsw_datum_ushort_vect_256("GEM_Engine_Map",val.map) ) );
  }

  template <>
  fsw_datum_inst<DFC_DB_stats>::fsw_datum_inst(const char* name, DFC_DB_stats& val)
    :fsw_datum(name),
     m_datum(val){
    //addChild( *(new fsw_datum_char_array_33_18("Stats_Titles",val.titles) ) );
  }

  template <>
  fsw_datum_inst<DFC_DB_prms>::fsw_datum_inst(const char* name, DFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<DFC_DB_gemV0>("DFC_GEM_map",val.gem) ) );
    addChild( *(new fsw_datum_inst<DFC_DB_stats>("DFC_Statistics",val.stats) ) );
  }

  template <>
  fsw_datum_inst<MFC_DB_calV0>::fsw_datum_inst(const char* name, MFC_DB_calV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_low" ,val.elo) ) );    
    addChild( *(new fsw_datum_uint("Energy_high",val.ehi) ) );    
    addChild( *(new fsw_datum_uint("nlayers"    ,val.nlayers) ) );    
  }

  template <>
  fsw_datum_inst<MFC_DB_tkrV0>::fsw_datum_inst(const char* name, MFC_DB_tkrV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("placeholder",val.placeholder) ) );   
  }

  template <>
  fsw_datum_inst<MFC_DB_prms>::fsw_datum_inst(const char* name, MFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<MFC_DB_calV0>("MFC_CAL",val.cal) ) );
    addChild( *(new fsw_datum_inst<MFC_DB_tkrV0>("MFC_TKR",val.tkr) ) );  
  }

  template <>
  fsw_datum_inst<HFC_DB_acdV0>::fsw_datum_inst(const char* name, HFC_DB_acdV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("placeholder",val.placeholder) ) ); 
  }

  template <>
  fsw_datum_inst<HFC_DB_calV0>::fsw_datum_inst(const char* name, HFC_DB_calV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Threshold" ,val.ethreshold) ) );    
    addChild( *(new fsw_datum_uint("n_Max_Logs"       ,val.nmaxlogs) ) );    
    addChild( *(new fsw_datum_uint("Tolerance"        ,val.tolerance) ) );        
  }

  template <>
  fsw_datum_inst<HFC_DB_tkrV0>::fsw_datum_inst(const char* name, HFC_DB_tkrV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("placeholder",val.placeholder) ) ); 
  }

  template <>
  fsw_datum_inst<HFC_DB_prms>::fsw_datum_inst(const char* name, HFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<HFC_DB_acdV0>("HFC_ACD",val.acd) ) );
    addChild( *(new fsw_datum_inst<HFC_DB_calV0>("HFC_CAL",val.cal) ) );  
    addChild( *(new fsw_datum_inst<HFC_DB_tkrV0>("HFC_TKR",val.tkr) ) );  
  }

  template <>
  fsw_datum_inst<GFC_DB_acdV0>::fsw_datum_inst(const char* name, GFC_DB_acdV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("TopSide_EnergyMax"       ,val.topSideEmax) ) );    
    addChild( *(new fsw_datum_uint("TopSideFilter_EnergyMax" ,val.topSideFilterEmax) ) );    
    addChild( *(new fsw_datum_uint("Splash_EnergyMax"        ,val.splashEmax) ) );        
    addChild( *(new fsw_datum_uint("Splash_Count"            ,val.splashCount) ) );    
    addChild( *(new fsw_datum_ptr("Splash_Map"               ,(void*)val.splashMap) ) );        
  }

  template <>
  fsw_datum_inst<GFC_DB_atfV0>::fsw_datum_inst(const char* name, GFC_DB_atfV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Max" ,val.emax) ) );    
    addChild( *(new fsw_datum_ptr("Shadowed_Map",(void*)val.shadowed) ) );            
  }

  template <>
  fsw_datum_inst<GFC_DB_zbottomV0>::fsw_datum_inst(const char* name, GFC_DB_zbottomV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Min" ,val.emin) ) );    
  }

  template <>
  fsw_datum_inst<GFC_DB_calV1>::fsw_datum_inst(const char* name, GFC_DB_calV1& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Pass"       ,val.epass) ) );    
    addChild( *(new fsw_datum_uint("Energy_Min"        ,val.emin) ) );    
    addChild( *(new fsw_datum_uint("Energy_Max"        ,val.emax) ) );        
    addChild( *(new fsw_datum_uint("Layer0_RatioLow"   ,val.layer0RatioLo) ) );    
    addChild( *(new fsw_datum_uint("Layer0_RatioHigh"  ,val.layer0RatioHi) ) );    
  }

  template <>
  fsw_datum_inst<GFC_DB_tkrV0>::fsw_datum_inst(const char* name, GFC_DB_tkrV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Row2_Energy_Max"       ,val.row2Emax) ) );    
    addChild( *(new fsw_datum_uint("Row1_Energy_Max"       ,val.row01Emax) ) );    
    addChild( *(new fsw_datum_uint("Top_Energy_Max"        ,val.topEmax) ) );        
    addChild( *(new fsw_datum_uint("Zero_Track_Energy_Min" ,val.zeroTkrEmin) ) );    
    addChild( *(new fsw_datum_uint("Two_Track_Energy_Max"  ,val.twoTkrEmax) ) );    
    addChild( *(new fsw_datum_uint("Skirt_Energy_Max"      ,val.skirtEmax) ) );    
  }

  template <>
  fsw_datum_inst<GFC_DB_prms>::fsw_datum_inst(const char* name, GFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<GFC_DB_acdV0>("GFC_ACD",val.acd) ) );
    addChild( *(new fsw_datum_inst<GFC_DB_atfV0>("GFC_ATF",val.atf) ) );
    addChild( *(new fsw_datum_inst<GFC_DB_zbottomV0>("GFC_ZBottom",val.zbottom) ) );
    addChild( *(new fsw_datum_inst<GFC_DB_calV1>("GFC_CAL",val.cal) ) );
    addChild( *(new fsw_datum_inst<GFC_DB_tkrV0>("GFC_TKR",val.tkr) ) );  
  }

  template <>
  fsw_datum_inst<GFC_DB_schema>::fsw_datum_inst(const char* name, GFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){   
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<GFC_DB_prms>("GFC_Params",val.prms) ) );    
  }

  template <>
  fsw_datum_inst<HFC_DB_schema>::fsw_datum_inst(const char* name, HFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<HFC_DB_prms>("HFC_Params",val.prms) ) );    
  }

  template <>
  fsw_datum_inst<MFC_DB_schema>::fsw_datum_inst(const char* name, MFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<MFC_DB_prms>("MFC_Params",val.prms) ) );    
  }

  template <>
  fsw_datum_inst<DFC_DB_schema>::fsw_datum_inst(const char* name, DFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<DFC_DB_prms>("DFC_Params",val.prms) ) );    
  }


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

  template <>
  void fsw_datum_inst<EFC_DB_prescalers>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' ' << get_name() << std::endl;
    for ( unsigned int i(0); i < 32; i++) {
      os << std::setw(depth) << ' ' 
	 << "  prescale"  << std::right
	 << std::setw(25-depth) << std::left << std::dec << i 
	 << ':' << m_datum.pv[i] << std::endl;
    }
  }
  
  template <>
  void fsw_datum_inst<DFC_DB_gemV0>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' ' << get_name() << std::endl;
    os << std::setw(depth) << ' ' 
       <<  "     :  0 :  1 :  2 :  3 :  4 :  5 :  6 :  7 :  8 :  9 :  a :  b :  c :  d :  e :  f" 
       << std::endl;
    unsigned int idx(0);
    for ( unsigned int i(0); i < 16; i++) {
      os << std::setw(depth) << ' ' << "  " << std::hex << std::setw(1) << std::setfill('0') << i << "   ";
      for ( unsigned int j(0); j < 16; j++, idx++) {
	os << std::hex << std::setw(4) << std::setfill('0') << m_datum.map[idx] << " ";
      }
      os << std::endl;
    }
  }

  template <>
  void fsw_datum_inst<DFC_DB_stats>::specializedPrint( int depth, std::ostream& os) const {
     os << std::setw(depth) << ' ' << get_name() << std::endl;
     for ( unsigned int i(0); i < 33; i++) {
       os << std::setw(depth) << ' ' << "  " << std::dec << std::setw(2) << i << ' ';
       const unsigned char* title = m_datum.titles[i];
       	if ( title == 0 ) {
	  os << std::endl;
	} else {
	  os << title << std::endl;
	}
     }
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

  template <>
  DOMElement* fsw_datum_inst<EFC_DB_prescalers>::specializedWriteToXml( DOMElement& parent ) const {
   DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
   if ( thisNode == 0 ) return 0;
   return thisNode;
  }

  template <>
  DOMElement* fsw_datum_inst<DFC_DB_gemV0>::specializedWriteToXml( DOMElement& parent ) const {
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
    if ( thisNode == 0 ) return 0;
    return thisNode;
  }
  
  template <>
  DOMElement* fsw_datum_inst<DFC_DB_stats>::specializedWriteToXml( DOMElement& parent ) const {
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
    if ( thisNode == 0 ) return 0;
    return thisNode;
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

  template <>
  bool fsw_datum_inst<EFC_DB_prescalers>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

  template <>
  bool fsw_datum_inst<DFC_DB_gemV0>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

  template <>
  bool fsw_datum_inst<DFC_DB_stats>::specializedReadFromXml( DOMElement& thisNode ) {
    return true;
  }

}
