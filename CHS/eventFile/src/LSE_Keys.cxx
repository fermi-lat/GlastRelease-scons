#include <stdio.h>

#include "eventFile/LSE_Keys.h"

namespace eventFile {

  void LSE_Keys::dump( const char* pre, const char* post ) const
  {
    printf( "%sLSE_Keys: LATC_master = 0x%08X%s", pre, LATC_master, post );
    printf( "%sLSE_Keys: LATC_ignore = 0x%08X%s", pre, LATC_ignore, post );
  }

  void LPA_Keys::dump( const char* pre, const char* post ) const
  {
    LSE_Keys::dump( pre, post );
    std::vector<unsigned>::const_iterator itr(CDM_keys.begin());
    int i=0;
    for ( ; itr!=CDM_keys.end(); ++itr, ++i ) {
      printf( "%sCMD_key[%d] = 0x%08X%s", pre, i, *itr, post );
    } 
  }

  void LCI_Keys::dump( const char* pre, const char* post ) const
  {
    LSE_Keys::dump( pre, post );
    printf( "%sLCI_script = 0x%08X%s", pre, LCI_script, post );
  }

}
