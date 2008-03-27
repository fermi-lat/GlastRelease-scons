//  -*- Mode: c++ -*-
//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Implementation of IO conversion for FSW CDM schema
//
//      This is the only file that needs to be updated to handle new FSW schema
//
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
#include "./fsw_iface.h"
#include "./fsw_enums.h"
#include "./fsw_cdb_inst.h"

// Other headers
#include "xmlBase/Dom.h"
#include "facilities/Util.h"

// FSW headers
// This is the only place these headers are included
#include "fsw/GFC_DB/GFC_DB_schema.h"
#include "fsw/XFC_DB/HFC_DB_schema.h"
#include "fsw/XFC_DB/MFC_DB_schema.h"
#include "fsw/XFC_DB/DFC_DB_schema.h"


namespace configData {

  // read a CDM and cast it to the correct type
  fsw_cdb* fsw_datum::read_cdb( const char* cdm_name, int option ) {
    
    // Bookkeeping about the cdm
    unsigned int key(0);
    short schemaID(-1);
    short versionID(-1);
    short instanceID(-1);

    // Load the CDM, set the bookkeeping data
    void* schema_ptr = const_cast<void*>(fsw_get_data( cdm_name, option,
						       &key, &schemaID, &versionID, &instanceID ));

    // Wrap and cast the loaded memory
    fsw_cdb* retVal(0);
    switch ( schemaID ) {
    case fsw_enums::GFC_SCHEMA_ID :
      retVal = new fsw_cdb_inst<GFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<GFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::DFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<DFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<DFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::MFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<MFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<MFC_DB_schema*>(schema_ptr));
      break;
    case fsw_enums::HFC_SCHEMA_ID:
      retVal = new fsw_cdb_inst<HFC_DB_schema>(cdm_name,key,schemaID,versionID,instanceID,
					       static_cast<HFC_DB_schema*>(schema_ptr));
      break;
    default:
      break;
    }
    return retVal;
  } 

  // Specializations    


  /**
   * Define the EFC_DB_prescalers structure
   *
   *   uint in      -> Input_Prescale
   *   uint out     -> Output_Prescale
   *   uint pv[32]  -> Prescale_XX
   *
   **/
  template <>
  fsw_datum_inst<EFC_DB_prescalers>::fsw_datum_inst(const char* name, EFC_DB_prescalers& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Input_Prescale",val.in) ) );
    addChild( *(new fsw_datum_uint("Output_Prescale",val.out) ) );
  }


  /**
   * Define the EFC_DB_sampler structure
   *
   *   mask               enabled      -> Enabled
   *   EFC_DB_prescalers  prescalers   -> Prescalers
   *
   **/
  template <>
  fsw_datum_inst<EFC_DB_sampler>::fsw_datum_inst(const char* name, EFC_DB_sampler& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_mask("Enabled",val.enabled) ) );
    addChild( *(new fsw_datum_inst<EFC_DB_prescalers>("Prescalers",val.prescalers) ) );
  }
    
  
  /**
   * Define the DFC_DB_gemV0 structure
   *
   *   uint               map[256]     -> GEM_Map_XX
   *
   **/
  template <>
  fsw_datum_inst<DFC_DB_gemV0>::fsw_datum_inst(const char* name, DFC_DB_gemV0& val)
    :fsw_datum(name),
     m_datum(val){
  }

  /**
   * Define the DFC_DB_stats structure
   *
   *   uchar              titles[33][18] -> Stats_Titles
   *
   **/
  template <>
  fsw_datum_inst<DFC_DB_stats>::fsw_datum_inst(const char* name, DFC_DB_stats& val)
    :fsw_datum(name),
     m_datum(val){
  }

  /**
   * Define the DFC_DB_prms structure
   *
   *   DFC_GEM_map         gem   -> DFC_GEM_map
   *   DFC_DB_stats        stats -> DFC_Statistics
   *
   **/
  template <>
  fsw_datum_inst<DFC_DB_prms>::fsw_datum_inst(const char* name, DFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<DFC_DB_gemV0>("DFC_GEM_map",val.gem) ) );
    addChild( *(new fsw_datum_inst<DFC_DB_stats>("DFC_Statistics",val.stats) ) );
  }

  /**
   * Define the MFC_DB_calV0 structure
   *
   *   uint elo     -> Energy_low
   *   uint ehi     -> Energy_high
   *   uint nlayers -> NLayers
   *
   **/
  template <>
  fsw_datum_inst<MFC_DB_calV0>::fsw_datum_inst(const char* name, MFC_DB_calV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_low" ,val.elo) ) );    
    addChild( *(new fsw_datum_uint("Energy_high",val.ehi) ) );    
    addChild( *(new fsw_datum_uint("N_Layers"   ,val.nlayers) ) );    
  }

  /**
   * Define the MFC_DB_tkrV0 structure
   *
   *   uint placeholder -> Placeholder
   *
   **/  
  template <>
  fsw_datum_inst<MFC_DB_tkrV0>::fsw_datum_inst(const char* name, MFC_DB_tkrV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Placeholder",val.placeholder) ) );   
  }

  /**
   * Define the MFC_DB_prms structure
   *
   *   MFC_DB_calV0 cal   -> MFC_CAL
   *   MFC_DB_tkrV0 tkr   -> MFC_TKR
   *
   **/  
  template <>
  fsw_datum_inst<MFC_DB_prms>::fsw_datum_inst(const char* name, MFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<MFC_DB_calV0>("MFC_CAL",val.cal) ) );
    addChild( *(new fsw_datum_inst<MFC_DB_tkrV0>("MFC_TKR",val.tkr) ) );  
  }

  /**
   * Define the HFC_DB_tkrV0 structure
   *
   *   uint placeholder -> Placeholder
   *
   **/  
  template <>
  fsw_datum_inst<HFC_DB_acdV0>::fsw_datum_inst(const char* name, HFC_DB_acdV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Placeholder",val.placeholder) ) ); 
  }

  /**
   * Define the HFC_DB_calV0 structure
   *
   *   uint ethreshold   -> Energy_Threshold
   *   uint nmaxlogs     -> N_Max_Logs
   *   uint tolerance    -> Tolerance
   *
   **/
  template <>
  fsw_datum_inst<HFC_DB_calV0>::fsw_datum_inst(const char* name, HFC_DB_calV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Threshold" ,val.ethreshold) ) );    
    addChild( *(new fsw_datum_uint("N_Max_Logs"       ,val.nmaxlogs) ) );    
    addChild( *(new fsw_datum_uint("Tolerance"        ,val.tolerance) ) );        
  }

  /**
   * Define the HFC_DB_tkrV0 structure
   *
   *   uint placeholder -> Placeholder
   *
   **/
  template <>
  fsw_datum_inst<HFC_DB_tkrV0>::fsw_datum_inst(const char* name, HFC_DB_tkrV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Placeholder",val.placeholder) ) ); 
  }

  /**
   * Define the HFC_DB_prms structure
   *
   *   HFC_DB_acdV0 acd   -> HFC_ACD
   *   HFC_DB_calV0 cal   -> HFC_CAL
   *   HFC_DB_tkrV0 tkr   -> HFC_TKR
   *
   **/  
  template <>
  fsw_datum_inst<HFC_DB_prms>::fsw_datum_inst(const char* name, HFC_DB_prms& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<HFC_DB_acdV0>("HFC_ACD",val.acd) ) );
    addChild( *(new fsw_datum_inst<HFC_DB_calV0>("HFC_CAL",val.cal) ) );  
    addChild( *(new fsw_datum_inst<HFC_DB_tkrV0>("HFC_TKR",val.tkr) ) );  
  }

  /**
   * Define the GFC_DB_acdV0 structure
   *
   *   uint  topSideEmax        -> TopSide_EnergyMax
   *   uint  topSideFilterEmax  -> TopSideFilter_EnergyMax
   *   uint  splashEmax         -> Splash_EnergyMax
   *   uint  splashCount        -> Splash_Count
   *   void* splashMap          -> Splash_Map
   *
   **/    
  template <>
  fsw_datum_inst<GFC_DB_acdV0>::fsw_datum_inst(const char* name, GFC_DB_acdV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("TopSide_Energy_Max"       ,val.topSideEmax) ) );    
    addChild( *(new fsw_datum_uint("TopSideFilter_Energy_Max" ,val.topSideFilterEmax) ) );    
    addChild( *(new fsw_datum_uint("Splash_Energy_Max"        ,val.splashEmax) ) );        
    addChild( *(new fsw_datum_uint("Splash_Count"             ,val.splashCount) ) );    
    addChild( *(new fsw_datum_ptr("Splash_Map"                ,(void*)val.splashMap) ) );        
  }

  /**
   * Define the GFC_DB_atfV0 structure
   *
   *   uint  emax               -> Energy_Max
   *   void* shadowed           -> Shadowed_Map
   *
   **/    
  template <>
  fsw_datum_inst<GFC_DB_atfV0>::fsw_datum_inst(const char* name, GFC_DB_atfV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Max" ,val.emax) ) );    
    addChild( *(new fsw_datum_ptr("Shadowed_Map",(void*)val.shadowed) ) );            
  }

  /**
   * Define the GFC_DB_zbottomV0 structure
   *
   *   uint  emin               -> Energy_Min
   *
   **/    
  template <>
  fsw_datum_inst<GFC_DB_zbottomV0>::fsw_datum_inst(const char* name, GFC_DB_zbottomV0& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Min" ,val.emin) ) );    
  }

  /**
   * Define the GFC_DB_calV1 structure
   *
   *   uint epass          -> Energy_Pass
   *   uint emin           -> Energy_Min
   *   uint emax           -> Energy_Max
   *   uint layer0RatioLo  -> Layer0_Ratio_Low
   *   uint layer0RatioHi  -> Layer0_Ratio_High
   *
   **/    
  template <>
  fsw_datum_inst<GFC_DB_calV1>::fsw_datum_inst(const char* name, GFC_DB_calV1& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_uint("Energy_Pass"       ,val.epass) ) );    
    addChild( *(new fsw_datum_uint("Energy_Min"        ,val.emin) ) );    
    addChild( *(new fsw_datum_uint("Energy_Max"        ,val.emax) ) );        
    addChild( *(new fsw_datum_uint("Layer0_Ratio_Low"  ,val.layer0RatioLo) ) );    
    addChild( *(new fsw_datum_uint("Layer0_Ratio_High" ,val.layer0RatioHi) ) );    
  }

  /**
   * Define the GFC_DB_tkrV0 structure
   *
   *   uint row2Emax          -> Row2_Energy_Max
   *   uint row01Emax         -> Row1_Energy_Max
   *   uint topEmax           -> Top_Energy_Max
   *   uint zeroTkrEmin       -> Zero_Track_Energy_Min
   *   uint twoTkrEmax        -> Two_Track_Energy_Max
   *   uint skirtEmax         -> Skirt_Energy_Max
   *
   **/    
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

  /**
   * Define the GFC_DB_prms structure
   *
   *   GFC_DB_acdV0     acd       -> GFC_ACD
   *   GFC_DB_atrV0     atf       -> GFC_ATF
   *   GFC_DB_zbottomV0 zbottom   -> GFC_ZBottom
   *   GFC_DB_calV1     cal       -> GFC_CAL
   *   GFC_DB_tkrV0     tkr       -> GFC_TKR
   *
   **/  
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

  /**
   * Define the GFC_DB_schema structure
   *
   *   EFC_DB_sampler   sampler   -> EFC_Sampler
   *   GFC_DB_prms      prms      -> GFC_Params
   *
   **/  
  template <>
  fsw_datum_inst<GFC_DB_schema>::fsw_datum_inst(const char* name, GFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){   
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<GFC_DB_prms>("GFC_Params",val.prms) ) );    
  }

  /**
   * Define the HFC_DB_schema structure
   *
   *   EFC_DB_sampler   sampler   -> EFC_Sampler
   *   HFC_DB_prms      prms      -> HFC_Params
   *
   **/  
  template <>
  fsw_datum_inst<HFC_DB_schema>::fsw_datum_inst(const char* name, HFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<HFC_DB_prms>("HFC_Params",val.prms) ) );    
  }

  /**
   * Define the MFC_DB_schema structure
   *
   *   EFC_DB_sampler   sampler   -> EFC_Sampler
   *   MFC_DB_prms      prms      -> MFC_Params
   *
   **/  
  template <>
  fsw_datum_inst<MFC_DB_schema>::fsw_datum_inst(const char* name, MFC_DB_schema& val)
    :fsw_datum(name),
     m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<MFC_DB_prms>("MFC_Params",val.prms) ) );    
  }

  /**
   * Define the DFC_DB_schema structure
   *
   *   EFC_DB_sampler   sampler   -> EFC_Sampler
   *   DFC_DB_prms      prms      -> DFC_Params
   *
   **/  
  template <>
  fsw_datum_inst<DFC_DB_schema>::fsw_datum_inst(const char* name, DFC_DB_schema& val)
    :fsw_datum(name),
    m_datum(val){
    addChild( *(new fsw_datum_inst<EFC_DB_sampler>("EFC_Sampler",val.sampler) ) );  
    addChild( *(new fsw_datum_inst<DFC_DB_prms>("DFC_Params",val.prms) ) );    
  }
  

  /**
   * Print the EFC_DB_prescalers structure
   * 32 lines of:
   *   
   *   prescaleXX         : value
   *
   **/  
  template <>
  void fsw_datum_inst<EFC_DB_prescalers>::specializedPrint( int depth, std::ostream& os) const {
    os << std::setw(depth) << ' ' << get_name() << std::endl;
    for ( unsigned int i(0); i < 32; i++) {
      os << std::setw(depth) << ' ' 
	 << "  prescale"  << std::right
	 << std::setw(25-depth) << std::left << std::dec << i 
	 << ": " << m_datum.pv[i] << std::endl;
    }
  }
  
  /**
   * Print the DFC_DB_gemV0 structure
   * 16 x 16 table of 16 bit masks write as 4 hex characters
   *
   **/  
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
  
  /**
   * Print the DFC_DB_stats structure
   * 33 lines of 
   *   0  Engine 0................
   *   etc.
   *
   **/  
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

  
  /**
   * Write the EFC_DB_prescalers structure to xml
   *   <Prescalers>
   *     <Prescale line="XX">value</Prescale>
   *
   * Note that:
   *    <Input_Prescale>value</Input_Prescale>
   *    <Output_Prescale>0</Output_Prescale>
   *
   * Are also written as they are child elements
   *
   **/  
  template <>
  DOMElement* fsw_datum_inst<EFC_DB_prescalers>::specializedWriteToXml( DOMElement& parent ) const {
    static const std::string line("line");
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
    if ( thisNode == 0 ) return 0;
    for ( unsigned int i(0); i < 32; i++) {
      DOMElement* cNode = XmlUtil::makeChildNodeWithContent(*thisNode, "Prescale", m_datum.pv[i]);
      xmlBase::Dom::addAttribute(cNode,line,i);
     if ( cNode == 0 ) return 0;
    }
   return thisNode;
  }

  /**
   * Write the DFC_DB_gemV0 structure to xml
   *   <DFC_GEM_map>
   *     <GEM_Map condition="0x00">0x4000</GEM_Map>
   *     
   **/
  template <>
  DOMElement* fsw_datum_inst<DFC_DB_gemV0>::specializedWriteToXml( DOMElement& parent ) const {
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
    if ( thisNode == 0 ) return 0;
    static const std::string condition("condition");
    for ( unsigned int i(0); i < 256; i++) {
      char tempName[10];
      sprintf(tempName,"0x%02x",i);
      char tempValue[8];
      sprintf(tempValue,"0x%04x",m_datum.map[i]);
      DOMElement* cNode = XmlUtil::makeChildNodeWithContent(*thisNode, "GEM_MAP", tempValue);
      if ( cNode == 0 ) return 0;
      xmlBase::Dom::addAttribute(cNode,condition,std::string(tempName));
    }
    return thisNode;
  }
  
  /**
   * Write the DFC_DB_stats structure to xml
   *   <DFC_DFC_Statistics>
   *     <Stats_Title line="0">GEM Analyzed.....</Stats_Title>
   *     
   **/
  template <>
  DOMElement* fsw_datum_inst<DFC_DB_stats>::specializedWriteToXml( DOMElement& parent ) const {
    static const std::string line("line");
    DOMElement* thisNode = XmlUtil::makeChildNode( parent, get_name().c_str() );
    if ( thisNode == 0 ) return 0;
    for ( unsigned int i(0); i < 33; i++) {
      char tempTitle[20];
      sprintf(tempTitle,"%s",m_datum.titles[i]);
      DOMElement* cNode = XmlUtil::makeChildNodeWithContent(*thisNode, "Stats_Title", tempTitle);
      if ( cNode == 0 ) return 0;
      xmlBase::Dom::addAttribute(cNode,line,i);
    }
    return thisNode;
  }

  
  // Read from XML
  template <>
  bool fsw_datum_inst<EFC_DB_prescalers>::specializedReadFromXml( DOMElement& thisNode ) {
    static std::string Prescale("Prescale");
    static std::string line("line");
    std::vector<DOMElement*> cVect;
    std::string txtVal;
    xmlBase::Dom::getChildrenByTagName(&thisNode,Prescale,cVect);
    for ( std::vector<DOMElement*>::const_iterator itr = cVect.begin(); 
	  itr != cVect.end(); itr++ ) {
      if ( ! XmlUtil::getAttribute(thisNode,"condition",txtVal ) ) return false;
      int i(-1);
      try { i = facilities::Util::stringToInt(txtVal); }
      catch ( std::exception& ) { return false; }
      if ( i < 0 || i > 31) return false;
      XmlUtil::getText(*(*itr),txtVal);
      if ( txtVal.empty() ) return false;
      try { m_datum.pv[i] = facilities::Util::stringToUnsigned(txtVal); }
      catch ( std::exception& ) { return false; }
    }
    return true;
  }

  template <>
  bool fsw_datum_inst<DFC_DB_gemV0>::specializedReadFromXml( DOMElement& thisNode ) {
    static const std::string GEM_MAP("GEM_MAP");
    std::vector<DOMElement*> cVect;
    
    std::string txtVal;
    xmlBase::Dom::getChildrenByTagName(&thisNode,GEM_MAP,cVect);
    for ( std::vector<DOMElement*>::const_iterator itr = cVect.begin(); 
	  itr != cVect.end(); itr++ ) {
      if ( ! XmlUtil::getAttribute(thisNode,"condition",txtVal ) ) return false;
      int i(-1);
      try { i = facilities::Util::stringToInt(txtVal); }
      catch ( std::exception& ) { return false; }
      if ( i < 0 || i > 255) return false;
      XmlUtil::getText(*(*itr),txtVal);
      if ( txtVal.empty() ) return false;
      try { m_datum.map[i] = facilities::Util::stringToUnsigned(txtVal); }
      catch ( std::exception& ) { return false; }
    }
    return true;
  }

  template <>
  bool fsw_datum_inst<DFC_DB_stats>::specializedReadFromXml( DOMElement& thisNode ) {
    static const std::string Stats_Title("Stats_Title");
    static const std::string line("line");
    std::vector<DOMElement*> cVect;
    std::string txtVal;
    xmlBase::Dom::getChildrenByTagName(&thisNode,Stats_Title,cVect);
    unsigned char** charPtr = (unsigned char**)(m_datum.titles);
    for ( std::vector<DOMElement*>::const_iterator itr = cVect.begin(); 
	  itr != cVect.end(); itr++ ) {
      if ( ! XmlUtil::getAttribute(thisNode,"condition",txtVal ) ) return false;
      int i(-1);
      try { i = facilities::Util::stringToInt(txtVal); }
      catch ( std::exception& ) { return false; }      
      if ( i < 0 || i > 32) return false;
      XmlUtil::getText(*(*itr),txtVal);
      if ( txtVal.empty() ) return false;
      for ( int j(0); j < 18; j++ ) {
	charPtr[i][j] = (unsigned char)(txtVal[j]);
      }
    }
    return true;
  }

}


