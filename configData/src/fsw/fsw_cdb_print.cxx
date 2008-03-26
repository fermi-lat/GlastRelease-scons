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
#include "./fsw_cdb_print.h"

#include <iomanip>


namespace configData {

  namespace fsw_cdb_print {
    
    int printMask( const char* name, unsigned int mask, std::ostream& os ) {
      os << name;
      unsigned int test(0x1);
      for ( unsigned int i(0); i < 32; i++) {
	if ( test & mask ) {
	  os << 'X';
	} else {
	  os << '.';
	}
	test *= 2;
      }
      os << std::endl;
      return 0;
    }

    int printValue( const char* name, unsigned int value, std::ostream& os ) {
      os << name << value << std::endl;
      return 0;
    }
    

    // Prescaler and sampler sub-structures   
    int print_efc_db_prescalers( const EFC_DB_prescalers& d, std::ostream& os){ 
      os << " EFC Prescalers" << std::endl;
      printValue( "  input              : " , d.in, os);
      printValue( "  output             : " , d.out, os);
      for ( unsigned int i(0); i < 32; i++) {
	os <<     "  prescale " << std::setw(2) << std::dec << i << "        : " << d.pv[i] << std::endl;
      }
      return 0;
    }

    int print_efc_db_sampler( const EFC_DB_sampler& d, std::ostream& os){ 
      os << " EFC" << std::endl;
      printMask(  "  enabled            : ", d.enabled, os);
      print_efc_db_prescalers( d.prescalers, os);
      return 0;
    }
   

    // Diagnositic filter sub-structures
    int print_dfc_db_gem_V0( const DFC_DB_gemV0& d, std::ostream& os){ 
      os << " DFG GEM V0" << std::endl;
      os <<   "     :  0 :  1 :  2 :  3 :  4 :  5 :  6 :  7 :  8 :  9 :  a :  b :  c :  d :  e :  f" 
	 << std::endl;
      unsigned int idx(0);
      for ( unsigned int i(0); i < 16; i++) {
	os << "  " << std::hex << std::setw(1) << std::setfill('0') << i << "   ";
	for ( unsigned int j(0); j < 16; j++, idx++) {
	  os << std::hex << std::setw(4) << std::setfill('0') << d.map[idx] << " ";
	}
	os << std::endl;
      }
      return 0;
    }
    
    int print_dfc_db_stats( const DFC_DB_stats& d, std::ostream& os){ 
      os << " DFC STATS" << std::endl;
      for ( unsigned int i(0); i < 33; i++) {
	const unsigned char* title = d.titles[i];
	os << "  " << std::dec << std::setw(2) << i << ' ';
	if ( title == 0 ) {
	  os << std::endl;
	} else {
	  os << title << std::endl;
	}
      }
      return 0;
    }

    int print_dfc_db_prms( const DFC_DB_prms& d, std::ostream& os){ 
      print_dfc_db_gem_V0(d.gem  ,os);
      print_dfc_db_stats(d.stats ,os);
      return 0;
    }


    // MIP filter sub-structures
    int print_mfc_db_cal_V0( const MFC_DB_calV0& d, std::ostream& os){ 
      os << " MFC CAL V0" << std::endl;
      printValue( "  elo                : " , d.elo, os);
      printValue( "  ehi                : " , d.ehi, os);
      printValue( "  nlayers            : " , d.nlayers, os);
      return 0;
    }
    
    int print_mfc_db_tkr_V0( const MFC_DB_tkrV0& d, std::ostream& os){ 
      os << " MFC TKR V0" << std::endl;
      printValue( "  placeholder        : " , d.placeholder, os);
      return 0;
    }

    int print_mfc_db_prms( const MFC_DB_prms& d, std::ostream& os){ 
      print_mfc_db_cal_V0(d.cal ,os);
      print_mfc_db_tkr_V0(d.tkr ,os);
      return 0;
    }
    

    // Heavy-ion filter sub-structures
    int print_hfc_db_acd_V0( const HFC_DB_acdV0& d, std::ostream& os){ 
      os << " HFC ACD V0" << std::endl;
      printValue( "  placeholder        : " , d.placeholder, os);     
      return 0;
    }

    int print_hfc_db_cal_V0( const HFC_DB_calV0& d, std::ostream& os){ 
      os << " HFC CAL V0" << std::endl;
      printValue( "  ethreshold         : " , d.ethreshold, os);
      printValue( "  nmaxlogs           : " , d.nmaxlogs, os);      
      printValue( "  tolerance          : " , d.tolerance, os);
      return 0;
    }
    
    int print_hfc_db_tkr_V0( const HFC_DB_tkrV0& d, std::ostream& os){ 
      os << " HFC TKR V0" << std::endl;
      printValue( "  placeholder        : " , d.placeholder, os);
      return 0;
    }

    int print_hfc_db_prms( const HFC_DB_prms& d, std::ostream& os){ 
      print_hfc_db_acd_V0(d.acd ,os);
      print_hfc_db_cal_V0(d.cal ,os);
      print_hfc_db_tkr_V0(d.tkr ,os);
      return 0;
    }
    

    // Gamma filter sub-structures
    int print_gfc_db_acd_V0( const GFC_DB_acdV0& d, std::ostream& os){ 
      os << " GFC ACD V0" << std::endl;
      printValue( "  topSideEmax        : " , d.topSideEmax, os);
      printValue( "  topSideFilterEmax  : " , d.topSideFilterEmax, os);
      printValue( "  splashEmax         : " , d.splashEmax, os);
      printValue( "  splashCount        : " , d.splashCount, os);
      os <<       "  splashMap          : " << d.splashMap << std::endl;
      return 0;
    }

    int print_gfc_db_atf_V0( const GFC_DB_atfV0& d, std::ostream& os){ 
      os << " GFC ATF V0" << std::endl;
      printValue( "  emax               : " , d.emax, os);
      os <<       "  shadow             : " << d.shadowed << std::endl;
      return 0;
    }

    int print_gfc_db_zbottom_V0( const GFC_DB_zbottomV0& d, std::ostream& os){ 
      os << " GFC ZBOTTOM V0 " << std::endl;
      printValue( "  emin               : " , d.emin, os);
      return 0;
    }

    //int print_gfc_db_cal_V0( const GFC_DB_calV0& d, std::ostream& os){ 
    //  return 0;
    //}

    int print_gfc_db_cal_V1( const GFC_DB_calV1& d, std::ostream& os){ 
      os << " GFC CAL V1 " << std::endl;
      printValue( "  epass              : " , d.epass, os);
      printValue( "  emin               : " , d.emin, os);
      printValue( "  emax               : " , d.emax, os);
      printValue( "  layer0RatioLo      : " , d.layer0RatioLo, os);
      printValue( "  layer0RatioHi      : " , d.layer0RatioHi, os);
      return 0;
    }
    
    int print_gfc_db_tkr_V0( const GFC_DB_tkrV0& d, std::ostream& os){ 
      os << " GFC TKR V0 " << std::endl;
      printValue( "  row2Emax           : " , d.row2Emax, os);
      printValue( "  row01Emax          : " , d.row01Emax, os);
      printValue( "  topEmax            : " , d.topEmax, os);
      printValue( "  zeroTkrEmin        : " , d.zeroTkrEmin, os);
      printValue( "  twoTkrEmax         : " , d.twoTkrEmax, os);
      printValue( "  skirtEmax          : " , d.skirtEmax, os);
      return 0;
    }

    int print_gfc_db_prms( const GFC_DB_prms& d, std::ostream& os){ 
      print_gfc_db_acd_V0(d.acd ,os);
      print_gfc_db_atf_V0(d.atf ,os);
      print_gfc_db_zbottom_V0(d.zbottom ,os);
      // print_gfc_db_cal_V0(d.cal ,os);
      print_gfc_db_cal_V1(d.cal ,os);
      print_gfc_db_tkr_V0(d.tkr ,os);
      return 0;
    }
    

    // Top level filter configuration structures

    // Gamma Filter
    int print_gfc( const GFC_DB_schema* schema, std::ostream& os ) {
      os << std::endl;
      os << "GFC " << std::endl;
      print_gfc_db_prms(schema->prms,os);
      print_efc_db_sampler(schema->sampler,os);
      return 0;
    }
    

    // Heavy-ion filter
    int print_hfc( const HFC_DB_schema* schema, std::ostream& os ) {
      os << std::endl;
      os << "HFC " << std::endl;
      print_hfc_db_prms(schema->prms,os);
      print_efc_db_sampler(schema->sampler,os);      
      return 0;
    }
    
    // MIP filter
    int print_mfc( const MFC_DB_schema* schema, std::ostream& os ) {
      os << std::endl;
      os << "MFC " << std::endl;
      print_mfc_db_prms(schema->prms,os);
      print_efc_db_sampler(schema->sampler,os);
      return 0;
    }
    
    // Diagnostic filter
    int print_dfc(const  DFC_DB_schema* schema, std::ostream& os ) {
      os << std::endl;
      os << "DFC " << std::endl;
      print_dfc_db_prms(schema->prms,os);
      print_efc_db_sampler(schema->sampler,os);
      return 0;
    }
  }
}
