#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>

#ifdef WIN32

#include <io.h>

#define ftruncate( a, b ) _chsize( (a), (b) )
#undef _fileno
#define fileno( a ) _fileno( (a) )

#else

#include <unistd.h>
#include <sys/types.h>

#endif

#include <sstream>
#include <stdexcept>

#include "eventFile/LSEWriter.h"

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Keys.h"

#include "facilities/Util.h"

namespace eventFile {

  LSEWriter::LSEWriter( const std::string& filename, unsigned runid )
    : m_name( filename ), m_hdr()
  {
    // stash the runid in the header
    m_hdr.m_runid = runid;

    // pick up the MOOT key/alias from the environment
    unsigned mootKey = 0xFFFFFFF0;
    const char* envbuf = getenv( "LSEWRITER_MOOTKEY" );
    if ( envbuf ) {
      mootKey = strtoul( envbuf, NULL, 0 );
    }
    const char* mootAlias = getenv( "LSEWRITER_MOOTALIAS" );
    if ( !mootAlias ) {
      mootAlias = "LSEWRITER_UNSET";
    }

    // set the MOOT key/alias values
    m_hdr.set_moot_key( mootKey );
    m_hdr.set_moot_alias( mootAlias );

    // expand any environment variables in the filename
    facilities::Util::expandEnvVar( &m_name );

    // open the specified file
    if ( ( m_FILE = fopen( m_name.c_str(), "wb" ) ) == NULL ) {
      std::ostringstream ess;
      ess << "LSEWriter::LSEWriter: error opening " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the file header
    writeHeader();
  }

  LSEWriter::~LSEWriter()
  {
    close();
  }

  void LSEWriter::close()
  {
    if ( m_FILE ) {
      writeHeader();
      if ( m_hdr.m_evtcnt == 0ULL ) {
	off_t zero(0);
	ftruncate( fileno( m_FILE ), zero );
      }
      fclose( m_FILE );
      m_FILE = NULL;
    }
  }

#ifdef _FILE_OFFSET_BITS
  off_t LSEWriter::tell()
  {
    return ftello( m_FILE );
  }

  void LSEWriter::writeHeader()
  {
    // seek to the beginning of the file
    off_t ofst = 0;
    fseeko( m_FILE, ofst, SEEK_SET );

    // write out the header data
    m_hdr.write( m_FILE );

    // return the location to the end of the file
    fseeko( m_FILE, ofst, SEEK_END );
  }
#else
  long LSEWriter::tell()
  {
    return ftell( m_FILE );
  }

  void LSEWriter::writeHeader()
  {
    // seek to the beginning of the file
    int ofst = 0;
    fseek( m_FILE, ofst, SEEK_SET );

    // write out the header data
    m_hdr.write( m_FILE );

    // return the location to the end of the file
    fseek( m_FILE, ofst, SEEK_END );
  }
#endif

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf )
  {
    size_t nitems(0);
    nitems = fwrite( reinterpret_cast<const void*>( &ctx ), sizeof( LSE_Context ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_File::write: error writing LSE_Context to " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    ebf.write( m_FILE );

    // capture header information
    if ( m_hdr.m_evtcnt == 0ULL ) {
      m_hdr.m_secs_beg = ctx.current.timeSecs;
      m_hdr.m_GEMseq_beg = ctx.scalers.sequence;
    }
    m_hdr.m_evtcnt++;
    m_hdr.m_secs_end = ctx.current.timeSecs;
    m_hdr.m_GEMseq_end = ctx.scalers.sequence;
  }
  
  void LSEWriter::write( int itype, const void* buf, size_t len )
  {
    // write the object type id to the file
    size_t nitems(0);
    nitems = fwrite( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LSE_Info typeid  to " << m_name;
      ess << " for type = " << itype << " size = " << len;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the object size to the file
    uint32_t const flen(static_cast<uint32_t>(len));
    nitems = fwrite( &flen, sizeof flen, 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LSE_Info size  to " << m_name;
      ess << " for type = " << itype << " size = " << len;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the object content to the file
    nitems = fwrite( buf, len, 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LSE_Info content to " << m_name;
      ess << " for type = " << itype << " size = " << len;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEWriter::write( const LPA_Keys& keys )
  {
    // write the object type id to the file
    size_t nitems(0);
    int itype = LSE_Keys::LPA;
    nitems = fwrite( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LPA_Keys typeid  to " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the two LSE keys and the LPA_db key
    unsigned ukeys[4];
    ukeys[0] = keys.LATC_master;
    ukeys[1] = keys.LATC_ignore;
    ukeys[2] = keys.SBS;
    ukeys[3] = keys.LPA_db;
    nitems = fwrite( ukeys, sizeof( ukeys ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LPA_Keys content to " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEWriter::write( const LCI_Keys& keys )
  {
    // write the object type id to the file
    size_t nitems(0);
    int itype = LSE_Keys::LCI;
    nitems = fwrite( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LCI_Keys typeid  to " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write the three unsigned keys to the file
    unsigned ukeys[3];
    ukeys[0] = keys.LATC_master;
    ukeys[1] = keys.LATC_ignore;
    ukeys[2] = keys.LCI_script;
    nitems = fwrite( ukeys, sizeof( ukeys ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEWriter::write: error writing LCI_Keys content to " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LPA_Info& info, const LPA_Keys& keys )
  {
    write( ctx, ebf );
    info.write( m_FILE );
    write( keys );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_ACD_Info& info, const LCI_Keys& keys )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_ACD;
    write( itype, &info, sizeof( info ) );
    write( keys );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_CAL_Info& info, const LCI_Keys& keys )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_CAL;
    write( itype, &info, sizeof( info ) );
    write( keys );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_TKR_Info& info, const LCI_Keys& keys )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_TKR;
    write( itype, &info, sizeof( info ) );
    write( keys );
  }

}
