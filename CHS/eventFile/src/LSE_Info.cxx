#include <stdio.h>

#include "eventFile/LSE_Info.h"

namespace eventFile {

  void LSE_Info::dump( const char* pre, const char* post ) const
  {
    printf( "%stimeTics = 0x%08X%s", pre, timeTics, post );
    timeHack.dump( pre, post );
  }

  void LPA_Info::dump() const
  {
    LSE_Info::dump( "LSE_Info: ", "\n" );
    printf( "LPA_Info: softwareKey = 0x%08X\n", softwareKey ); 
    printf( "LPA_Info: hardwareKey = 0x%08X\n", hardwareKey ); 
  }

  void LCI_Info::dump( const char* pre, const char* post ) const
  {
    printf( "%sautoRange        = %s%s", pre, (autoRange) ? "True" : "False", post );
    printf( "%szeroSupression   = %s%s", pre, (zeroSupression) ? "True" : "False", post );
    printf( "%speriodicPrescale = %10d%s", pre, periodicPrescale, post );
    printf( "%ssoftwareKey      = 0x%08X%s", pre, softwareKey, post );
    printf( "%swriteCfg         = 0x%08X%s", pre, writeCfg, post );
    printf( "%sreadCfg          = 0x%08X%s", pre, readCfg, post );
  }

  void LCI_Channel::dump( const char* pre, const char* post ) const
  {
    printf( "%ssingle = %6d%s", pre, single, post );
    printf( "%sall    = %s%s", pre, (all) ? "True" : "False", post );
    printf( "%slatc   = %s%s", pre, (latc) ? "True" : "False", post );
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
    printf( "%she = 0x%04X%s", pre, he, post );
  }

  void LCI_ACD_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_ACD_Info: biasDac   = 0x%04X\n", biasDac );
    printf( "LCI_ACD_Info: injected  = 0x%04X\n", injected );
    printf( "LCI_ACD_Info: threshold = 0x%04X\n", threshold );
    printf( "LCI_ACD_Info: holdDelay = 0x%04X\n", holdDelay );
    trigger.dump( "LCI_AcdTrigger: ", "\n" );
    channel.dump( "LCI_Channel: ", "\n" );
  }

  void LCI_CAL_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_CAL_Info: uld       = 0x%04X\n", uld );
    printf( "LCI_CAL_Info: injected  = 0x%04X\n", injected );
    printf( "LCI_CAL_Info: threshold = 0x%04X\n", threshold );
    printf( "LCI_CAL_Info: delay     = 0x%04X\n", delay );
    trigger.dump( "LCI_CalTrigger: ", "\n" );
    channel.dump( "LCI_Channel: ", "\n" );
  }

  void LCI_TKR_Info::dump() const
  {
    LCI_Info::dump( "LCI_Info: ", "\n" );
    printf( "LCI_TKR_Info: injected  = 0x%04X\n", injected );
    printf( "LCI_TKR_Info: threshold = 0x%04X\n", threshold );
    printf( "LCI_TKR_Info: delay     = 0x%04X\n", delay );
    channel.dump( "LCI_Channel: ", "\n" );
  }
}
