// -*- Mode: c++ -*-
#ifndef fsw_cdb_print_HH
#define fsw_cdb_print_HH
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
/** @namespace fsw_cdb_print
    @author E. Charles

    Print functions to convert FSW CDM to xml
*/

// c++/stl include
#include <iostream>

// FSW includes
#include "GFC_DB/GFC_DB_schema.h"
#include "XFC_DB/HFC_DB_schema.h"
#include "XFC_DB/MFC_DB_schema.h"
#include "XFC_DB/DFC_DB_schema.h"


namespace configData {
  
  namespace fsw_cdb_print {

    int printMask( const char* name, unsigned int mask, std::ostream& os = std::cout );

    int printValue( const char* name, unsigned int value, std::ostream& os = std::cout );
    

    // Prescaler and sampler sub-structures   
    int print_efc_db_prescalers( const EFC_DB_prescalers& d, std::ostream& os = std::cout );

    int print_efc_db_sampler( const EFC_DB_sampler& d, std::ostream& os = std::cout );
   

    // Diagnositic filter sub-structures
    int print_dfc_db_gem_V0( const DFC_DB_gemV0& d, std::ostream& os = std::cout );
    
    int print_dfc_db_stats( const DFC_DB_stats& d, std::ostream& os = std::cout );

    int print_dfc_db_prms( const DFC_DB_prms& d, std::ostream& os = std::cout );


    // MIP filter sub-structures
    int print_mfc_db_cal_V0( const MFC_DB_calV0& d, std::ostream& os = std::cout );
    
    int print_mfc_db_tkr_V0( const MFC_DB_tkrV0& d, std::ostream& os = std::cout );

    int print_mfc_db_prms( const MFC_DB_prms& d, std::ostream& os = std::cout );
    

    // Heavy-ion filter sub-structures
    int print_hfc_db_acd_V0( const HFC_DB_acdV0& d, std::ostream& os = std::cout );

    int print_hfc_db_cal_V0( const HFC_DB_calV0& d, std::ostream& os = std::cout );
    
    int print_hfc_db_tkr_V0( const HFC_DB_tkrV0& d, std::ostream& os = std::cout );

    int print_hfc_db_prms( const HFC_DB_prms& d, std::ostream& os = std::cout );
    

    // Gamma filter sub-structures
    int print_gfc_db_acd_V0( const GFC_DB_acdV0& d, std::ostream& os = std::cout );

    int print_gfc_db_atf_V0( const GFC_DB_atfV0& d, std::ostream& os = std::cout );

    int print_gfc_db_zbottom_V0( const GFC_DB_zbottomV0& d, std::ostream& os = std::cout );

    //int print_gfc_db_cal_V0( const GFC_DB_calV0& d, std::ostream& os = std::cout );

    int print_gfc_db_cal_V1( const GFC_DB_calV1& d, std::ostream& os = std::cout );
    
    int print_gfc_db_tkr_V0( const GFC_DB_tkrV0& d, std::ostream& os = std::cout );

    int print_gfc_db_prms( const GFC_DB_prms& d, std::ostream& os = std::cout );
    

    // Top level filter configuration structures
    int print_gfc( const GFC_DB_schema* schema, std::ostream& os = std::cout );
    
    int print_hfc( const HFC_DB_schema* schema, std::ostream& os = std::cout );
    
    int print_mfc( const MFC_DB_schema* schema, std::ostream& os = std::cout );
    
    int print_dfc( const DFC_DB_schema* schema, std::ostream& os = std::cout );

  }

}


#endif
