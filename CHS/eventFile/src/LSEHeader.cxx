#include <errno.h>
#include <string.h>

#include <sstream>
#include <stdexcept>

#include "eventFile/LSEHeader.h"

namespace eventFile {

  LSEHeader::LSEHeader() :
    m_runid(0), m_secs_beg(0), m_secs_end(0), m_evtcnt(0), m_GEMseq_beg(0), m_GEMseq_end(0)
  {
    memset( m_src_apids,  0, LSEHEADER_MAX_APIDS*sizeof(unsigned) );
    memset( m_src_seqerr, 0, LSEHEADER_MAX_APIDS*sizeof(unsigned) );
    memset( m_dfi_apids,  0, LSEHEADER_MAX_APIDS*sizeof(unsigned) );
    memset( m_dfi_dfierr, 0, LSEHEADER_MAX_APIDS*sizeof(unsigned) );
  }

  void LSEHeader::read( FILE* fp )
  {
    // check for the correct file marker value
    unsigned marker( 0 );
    fread( &marker, sizeof(unsigned), 1, fp );
    if ( marker != 0xFAF32000 ) {
      std::ostringstream ess;
      ess << "LSEHeader::read: invalid header";
      throw std::runtime_error( ess.str() );
    }

    // read in the header data
    size_t nitems(0);
    nitems = fread( this, sizeof( LSEHeader ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEHeader::read: error reading header, ";
      ess << "(" << errno << ":'" << strerror( errno ) << "')";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEHeader::write( FILE* fp )
  {
    // write out the header-marker value 0xFAF32000
    unsigned marker = 0xFAF32000;
    fwrite( &marker, sizeof(unsigned), 1, fp );

    // write out the header data
    fwrite( this, sizeof( LSEHeader ), 1, fp );
  }
}
