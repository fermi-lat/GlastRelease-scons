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

}
