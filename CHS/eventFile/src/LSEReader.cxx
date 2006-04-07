#include <errno.h>

#include <sstream>
#include <stdexcept>

#include "eventFile/LSEReader.h"

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/EBF_Data.h"

#include "facilities/Util.h"

namespace eventFile {

  LSEReader::LSEReader( const std::string& filename )
    : m_name( filename ), m_runid( 0 ),
      m_evtcnt( 0 ), m_GEMseq_beg( 0 ), m_GEMseq_end( 0 )
  {
    // expand any environment variables in the filename
    facilities::Util::expandEnvVar( &m_name );

    // open the specified file
    if ( ( m_FILE = fopen( m_name.c_str(), "rb" ) ) == NULL ) {
      std::ostringstream ess;
      ess << "LSEReader::LSEReader: error opening " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // read in the file header
    readHeader();
  }

  LSEReader::~LSEReader()
  {
    close();
  }

  void LSEReader::close()
  {
    if ( m_FILE ) {
      fclose( m_FILE );
      m_FILE = NULL;
    }
  }

  void LSEReader::readHeader()
  {
    // seek to the beginning of the file
    off_t ofst = 0;
    fseeko( m_FILE, ofst, SEEK_SET ); 

    // check for the correct file marker value
    unsigned marker(0);
    fread( &marker, sizeof(unsigned), 1, m_FILE );
    if ( marker != 0xFAF32000 ) {
      std::ostringstream ess;
      ess << "LSEReader::readHeader: invalid header in " << m_name;
      ess << "; not an LSE file";
      throw std::runtime_error( ess.str() );
    }

    // read in the run ID
    fread( &m_runid, sizeof(unsigned), 1, m_FILE );

    // read in the event statistics
    fread( &m_evtcnt,     sizeof(unsigned long long), 1, m_FILE );
    fread( &m_GEMseq_beg, sizeof(unsigned long long), 1, m_FILE );
    fread( &m_GEMseq_end, sizeof(unsigned long long), 1, m_FILE );
  }

  bool LSEReader::read( LSE_Context& ctx, EBF_Data& ebf )
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
	ess << "LSEReader::read: error reading LSE_Context from " << m_name;
	ess << " (" << errno << "=" << strerror( errno ) << ")";
	throw std::runtime_error( ess.str() );
      }
    }

    // populate the supplied context object
    ctx = *( reinterpret_cast< LSE_Context* >( buf ) );

    // read in the EBF data
    ebf.read( m_FILE );

    return true;
  }

  void LSEReader::read( int& itype, unsigned char* buf, size_t& len )
  {
    // read in the LSE_Info typeid
    size_t nitems(0);
    nitems = fread( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info typeid from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // read in the LSE_Info size
    nitems = fread( &len, sizeof( size_t ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info size from " << m_name;
      ess << " for typeid = " << itype;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // read in the LSE_Info content
    nitems = fread( buf, len, 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info content from " << m_name;
      ess << " for typeid = " << itype << " size = " << len;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  bool LSEReader::read( LSE_Context& ctx, EBF_Data& ebf, LSE_Info::InfoType& infotype,
			LPA_Info& pinfo, LCI_ACD_Info& ainfo, LCI_CAL_Info& cinfo, LCI_TKR_Info& tinfo )
  {
    // read the context and EBF 
    if ( !read( ctx, ebf ) ) {
      return false;
    }

    // read the LSE_Info object content
    int itype(0);
    unsigned char buf[ 128 * 1024 ];
    size_t len(0);
    read( itype, &buf[0], len );

    // assign to the proper type of object
    infotype = (LSE_Info::InfoType) itype;
    switch ( infotype ) {
    case LSE_Info::LPA:
      pinfo = *reinterpret_cast<LPA_Info*>( buf );
      break;
    case LSE_Info::LCI_ACD:
      ainfo = *reinterpret_cast<LCI_ACD_Info*>( buf );
      break;
    case LSE_Info::LCI_CAL:
      cinfo = *reinterpret_cast<LCI_CAL_Info*>( buf );
      break;
    case LSE_Info::LCI_TKR:
      tinfo = *reinterpret_cast<LCI_TKR_Info*>( buf );
      break;
    }

    return true;
  }

}
