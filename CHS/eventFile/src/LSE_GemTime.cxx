#include <stdio.h>

#include "eventFile/LSE_GemTime.h"

namespace eventFile {

  void LSE_GemTime::dump( const char* pre, const char* post ) const
  {
    printf( "%stics  = 0x%08X (%d)%s", pre, tics, tics, post );
    printf( "%shacks = 0x%08X (%d)%s", pre, hacks, hacks, post );
  }

}
