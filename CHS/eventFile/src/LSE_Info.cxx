#include <errno.h>
#include <stdexcept>
#include <sstream>

#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"

namespace eventFile {

  void LSE_Info::dump( const char* pre, const char* post ) const
  {
    printf( "%scompressionLevel = 0x%08x%s", pre, compressionLevel, post );
    printf( "%scompressedSize   = 0x%08x (%d)%s", pre, compressedSize, compressedSize, post );
    printf( "%stimeTics         = 0x%08X%s", pre, timeTics, post );
    timeHack.dump( pre, post );
  }

  void LPA_Info::dump() const
  {
    LSE_Info::dump( "LSE_Info: ", "\n" );
    printf( "LPA_Info: softwareKey = 0x%08X\n", softwareKey ); 
    printf( "LPA_Info: hardwareKey = 0x%08X\n", hardwareKey ); 
    printf( "LPA_Info: lpaDbKey    = 0x%08X\n", lpaDbKey );
    std::vector< LPA_Handler >::const_iterator itr( handlers.begin() );
    int ihnd=0;
    for ( ; itr != handlers.end(); ++itr, ++ihnd ) {
      char pre[256];
      sprintf( pre, "LPA_Info: LPA_Handler[%X]: ", ihnd );
      itr->dump( pre, "\n" );
    }
  }

  void LPA_Info::write( FILE* fp ) const
  {
    // write the object type id to the file
    size_t nitems(0);
    int itype = LSE_Info::LPA;
    nitems = fwrite( &itype, sizeof( int ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::write: error writing LPA_Info typeid ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the "fixed" part of the structure to the file
    size_t fixedsize = sizeof( LPA_Info ) - sizeof( std::vector< LPA_Handler > );
    nitems = fwrite( this, fixedsize, 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::write: error writing fixed LPA_Info content ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the number of LPA_Handler instances to the file
    unsigned nhandlers = handlers.size();
    nitems = fwrite( &nhandlers, sizeof( unsigned ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LPA_Info nhandlers  ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // if there are no handlers, then we're done
    if ( nhandlers == 0 ) {
      return;
    }

    // write each LPA_Handler to the file
    nitems = fwrite( &(handlers[0]), nhandlers * sizeof( LPA_Handler ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::write: error writing LPA_Handler content ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LPA_Info::read( FILE* fp )
  {
    // read the fixed-size LPA_Info data from the file
    size_t fixedsize = sizeof( LPA_Info ) - sizeof( std::vector< LPA_Handler > );
    int nitems(0);
    nitems = fread( this, fixedsize, 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::read: error reading fixed LPA_Info content ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    
    // read the number of LPA_Handler instances
    unsigned nhandlers;
    nitems = fread( &nhandlers, sizeof( nhandlers ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::read: error reading handler count ";
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // if there are no handlers, then we're done
    if ( nhandlers == 0 ) {
      return;
    }

    // read the handler instances from the file
    handlers.resize( nhandlers );
    nitems = fread( &(handlers[0]), nhandlers * sizeof( LPA_Handler ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_Info::read: error reading handlers block ";
      ess << "nitems = " << nitems << " nhandlers = " << nhandlers;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LCI_Info::dump( const char* pre, const char* post ) const
  {
    LSE_Info::dump( "LSE_Info: ", "\n" );
    printf( "%sautoRange        = %s%s", pre, (autoRange) ? "True" : "False", post );
    printf( "%szeroSupression   = %s%s", pre, (zeroSupression) ? "True" : "False", post );
    printf( "%speriodicPrescale = %10d%s", pre, periodicPrescale, post );
    printf( "%ssoftwareKey      = 0x%08X%s", pre, softwareKey, post );
    printf( "%swriteCfg         = 0x%08X%s", pre, writeCfg, post );
    printf( "%sreadCfg          = 0x%08X%s", pre, readCfg, post );
  }

  void LCI_Channel::dump( const char* pre, const char* post ) const
  {
    printf( "%snumChan = 0x%04X%s", pre, numChan, post );
    printf( "%ssingle  = %s%s", pre, (single) ? "True" : "False", post );
    printf( "%sall     = %s%s", pre, (all) ? "True" : "False", post );
    printf( "%slatc    = %s%s", pre, (latc) ? "True" : "False", post );
    printf( "%sperFE   = %s%s", pre, (perFE) ? "True" : "False", post );
  }

  void LCI_AcdTrigger::dump( const char* pre, const char* post ) const
  {
    printf( "%sveto        = 0x%04X%s", pre, veto, post );
    printf( "%svetoVernier = 0x%04X%s", pre, vetoVernier, post );
    printf( "%shld         = 0x%04X%s", pre, hld, post );
  }

  void LCI_CalTrigger::dump( const char* pre, const char* post ) const
  {
    printf( "%sle = 0x%04X%s", pre, le, post );
    printf( "%slowTrgEna = %s%s", pre, (lowTrgEna == 0xffff) ? "Undefined" : (lowTrgEna) ?  "True" : "False", post );
    printf( "%she = 0x%04X%s", pre, he, post );
    printf( "%shighTrgEna = %s%s", pre, (highTrgEna == 0xffff) ? "Undefined" : (highTrgEna) ?  "True" : "False", post );
  }

  void LCI_ACD_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_ACD_Info: biasDac   = 0x%04X\n", biasDac );
    printf( "LCI_ACD_Info: injected  = 0x%04X\n", injected );
    printf( "LCI_ACD_Info: threshold = 0x%04X\n", threshold );
    printf( "LCI_ACD_Info: holdDelay = 0x%04X\n", holdDelay );
    printf( "LCI_ACD_Info: hitmapDelay = 0x%04X\n", hitmapDelay );
    printf( "LCI_ACD_Info: range = 0x%04X\n", range );
    trigger.dump( "LCI_AcdTrigger: ", "\n" );
    channel.dump( "LCI_Channel: ", "\n" );
  }

  void LCI_CAL_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_CAL_Info: uld        = 0x%04X\n", uld );
    printf( "LCI_CAL_Info: injected   = 0x%04X\n", injected );
    printf( "LCI_CAL_Info: threshold  = 0x%04X\n", threshold );
    printf( "LCI_CAL_Info: delay      = 0x%04X\n", delay );
    printf( "LCI_CAL_Info: firstRange = 0x%04X\n", firstRange );
    printf( "LCI_CAL_Info: calibGain  = %s\n",  (calibGain == 0xffff) ? "Undefined" : (calibGain) ?  "True" : "False");
    printf( "LCI_CAL_Info: highCalEna = %s\n",  (highCalEna == 0xffff) ? "Undefined" : (highCalEna) ?  "True" : "False");
    printf( "LCI_CAL_Info: highRngEna = %s\n",  (highRngEna == 0xffff) ? "Undefined" : (highRngEna) ?  "True" : "False");
    printf( "LCI_CAL_Info: highGain   = 0x%04X\n", highGain );
    printf( "LCI_CAL_Info: lowCalEna  = %s\n",  (lowCalEna == 0xffff) ? "Undefined" : (lowCalEna) ?  "True" : "False");
    printf( "LCI_CAL_Info: lowRngEna  = %s\n",  (lowRngEna == 0xffff) ? "Undefined" : (lowRngEna) ?  "True" : "False");
    printf( "LCI_CAL_Info: lowGain    = 0x%04X\n", lowGain );
    trigger.dump( "LCI_CalTrigger: ", "\n" );
    channel.dump( "LCI_Channel: ", "\n" );
  }

  void LCI_TKR_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_TKR_Info: dacRange  = 0x%04X\n", dacRange );
    printf( "LCI_TKR_Info: injected  = 0x%04X\n", injected );
    printf( "LCI_TKR_Info: threshold = 0x%04X\n", threshold );
    printf( "LCI_TKR_Info: delay     = 0x%04X\n", delay );
    printf( "LCI_TKR_Info: splitLow  = 0x%04X\n", splitLow );
    printf( "LCI_TKR_Info: splitHigh = 0x%04X\n", splitHigh );
    channel.dump( "LCI_Channel: ", "\n" );
  }
}
