#include <stdio.h>

#include "eventFile/LSE_Keys.h"

namespace eventFile {

  void LSE_Keys::dump( const char* pre, const char* post ) const
  {
    printf( "%sLATC_master = 0x%08X%s", pre, LATC_master, post );
    printf( "%sLATC_ignore = 0x%08X%s", pre, LATC_ignore, post );
  }

  void LPA_Keys::dump( const char* pre, const char* post ) const
  {
    LSE_Keys::dump( "LSE_Keys", "\n" );
    std::vector<unsigned>::const_iterator itr(LPA_DB.begin());
    int i=0;
    for ( ; itr!=LPA_DB.end(); ++itr, ++i ) {
      printf( "%sLPA_DB[%d] = 0x%08X%s", pre, i, *itr, post );
    } 
  }

  void LCI_Keys::dump( const char* pre, const char* post ) const
  {
    LSE_Keys::dump( "LSE_Keys", "\n" );
    printf( "%sLCI_script = 0x%08X%s", pre, LCI_script, post );
  }

}
