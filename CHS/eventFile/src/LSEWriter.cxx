#include <errno.h>

#include <sstream>
#include <stdexcept>

#include "eventFile/LSEWriter.h"

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/EBF_Data.h"

#include "facilities/Util.h"

namespace eventFile {

  LSEWriter::LSEWriter( const std::string& filename, unsigned runid )
    : m_name( filename ), m_hdr()
  {
    // stash the runid in the header
    m_hdr.m_runid = runid;

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
      fclose( m_FILE );
      m_FILE = NULL;
    }
  }

#ifdef _FILE_OFFSET_BITS
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
#warning "64-bit file operations unavailable, file size limited to 2GB"
  void LSEWriter::writeHeader()
  {
    // seek to the beginning of the file
    long ofst = 0;
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
    nitems = fwrite( &len, sizeof( size_t ), 1, m_FILE );
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

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LPA_Info& info )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LPA;
    write( itype, &info, sizeof( info ) );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_ACD_Info& info )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_ACD;
    write( itype, &info, sizeof( info ) );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_CAL_Info& info )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_CAL;
    write( itype, &info, sizeof( info ) );
  }

  void LSEWriter::write( const LSE_Context& ctx, const EBF_Data& ebf, const LCI_TKR_Info& info )
  {
    write( ctx, ebf );
    int itype = LSE_Info::LCI_TKR;
    write( itype, &info, sizeof( info ) );
  }

}
