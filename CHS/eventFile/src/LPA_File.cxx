#include <errno.h>

#include <sstream>
#include <stdexcept>

#include "eventFile/LPA_File.h"

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/EBF_Data.h"

#include "facilities/Util.h"

static const char* ModeStr[] = { "rb", "wb" };

namespace eventFile {

  LPA_File::LPA_File( const std::string& filename, LPA_File::Mode openmode, unsigned runid )
    : m_name( filename ), m_mode( ModeStr[openmode] ), m_runid( runid ),
      m_evtcnt( 0 ), m_GEMseq_beg( 0 ), m_GEMseq_end( 0 )
  {
    // expand any environment variables in the filename
    facilities::Util::expandEnvVar( &m_name );

    // open the specified file
    if ( ( m_FILE = fopen( m_name.c_str(), m_mode.c_str() ) ) == NULL ) {
      std::ostringstream ess;
      ess << "LPA_File::LPA_File: error opening " << m_name;
      ess << " with mode '" << m_mode << "' ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // if we're reading, read the runid; if writing, write the runid
    switch ( openmode ) {
    case Read:
      fread( &m_runid, sizeof( unsigned ), 1, m_FILE );
      break;
    case Write:
      fwrite( &m_runid, sizeof( unsigned ), 1, m_FILE );
      break;
    }
  }

  LPA_File::~LPA_File()
  {
    close();
  }

  void LPA_File::close()
  {
    if ( m_FILE ) {
      fclose( m_FILE );
      m_FILE = NULL;
    }
  }

  void LPA_File::write( const LSE_Context& ctx, const LPA_Info& info, const EBF_Data& ebf )
  {
    size_t nitems(0);
    nitems = fwrite( reinterpret_cast<const void*>( &ctx ), sizeof( LSE_Context ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_File::write: error writing LSE_Context to " << m_name;
      ess << " with mode '" << m_mode << "' ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    nitems = fwrite( reinterpret_cast<const void*>( &info ), sizeof( LPA_Info ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_File::write: error writing LPA_Info to " << m_name;
      ess << " with mode '" << m_mode << "' ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    ebf.write( m_FILE );
  }

  bool LPA_File::read( LSE_Context& ctx, LPA_Info& info, EBF_Data& ebf )
  {
    unsigned char buf[ 128 * 1024 ];
    size_t nitems(0);

    // see if we're at the end of the file
    if ( feof( m_FILE ) ) return false;

    // read the context data as a bag-o-bytes
    nitems = fread( buf, sizeof( LSE_Context ), 1, m_FILE );
    if ( nitems != 1 ) {
      if ( feof( m_FILE ) ) {
	return false;
      } else {
	std::ostringstream ess;
	ess << "LPA_File::read: error reading LSE_Context from " << m_name;
	ess << " with mode '" << m_mode << "' ";
	ess << "(" << errno << "=" << strerror( errno ) << ")";
	throw std::runtime_error( ess.str() );
      }
    }

    // populate the supplied context object
    ctx = *( reinterpret_cast< LSE_Context* >( buf ) );

    // read the meta-information as a bag-o-bytes
    nitems = fread( buf, sizeof( LPA_Info ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LPA_File::read: error reading LPA_Info from " << m_name;
      ess << " with mode '" << m_mode << "' ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // populate the supplied meta-information object
    info = *( reinterpret_cast< LPA_Info* >( buf ) );

    // read in the EBF data
    ebf.read( m_FILE );

    return true;
  }

}
