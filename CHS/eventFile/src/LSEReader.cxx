#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <cstring>

#include <sstream>
#include <stdexcept>

#ifdef HAVE_FACILITIES
#include "facilities/Util.h"
#endif

#include "eventFile/LSEReader.h"
#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"
#include "eventFile/EBF_Data.h"

namespace eventFile {

  LSEReader::LSEReader( const std::string& filename )
    : m_name( filename ), m_hdr()
  {
#ifdef HAVE_FACILITIES
    // expand any environment variables in the filename
    facilities::Util::expandEnvVar( &m_name );
#endif

#ifndef _FILE_OFFSET_BITS
    // on windows, check to see if the file is >2GB in size and
    // throw an exception if so.
    struct _stati64 stbuf;
    if ( _stati64( m_name.c_str(), &stbuf ) == 0 ) {
      if ( stbuf.st_size > 0x7FFFFFFFLL ) {
	std::ostringstream ess;
	ess << "LSEReader::LSEReader: " << m_name;
	ess << " too large to open on Windows";
	throw std::runtime_error( ess.str() );
      }
    }
#endif

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

#ifdef _FILE_OFFSET_BITS
  void LSEReader::readHeader()
  {
    // seek to the beginning of the file
    off_t ofst = 0;
    fseeko( m_FILE, ofst, SEEK_SET ); 

    // read in the header data
    m_hdr.read( m_FILE );
  }

  int LSEReader::seek( off_t ofst )
  {
    return fseeko( m_FILE, ofst, SEEK_SET );
  }
#else
  void LSEReader::readHeader()
  {
    // seek to the beginning of the file
    int ofst = 0;
    fseek( m_FILE, ofst, SEEK_SET ); 

    // read in the header data
    m_hdr.read( m_FILE );
  }

  int LSEReader::seek( int ofst )
  {
    return fseek( m_FILE, ofst, SEEK_SET );
  }
#endif

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

  void LSEReader::read( unsigned char* buf, size_t& len )
  {
    // read in the LSE_Info size
    int nitems(0);
    uint32_t flen(0);
    nitems = fread( &flen, sizeof flen, 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info size from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    len = flen;

    // read in the LSE_Info content
    nitems = fread( buf, len, 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info content from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEReader::read( LSE_Keys& keys )
  {
    // every keys object has LATC_master and LATC_ignore
    size_t nitems(0);
    nitems = fread( &keys.LATC_master, sizeof( unsigned ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Keys LATC_master from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    nitems = fread( &keys.LATC_ignore, sizeof( unsigned ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Keys LATC_ignore from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEReader::read( LPA_Keys& pakeys )
  {
    // read the common part of the struct
    LSE_Keys& ekeys = pakeys;
    read( ekeys );

    // read the SBS value
    size_t nitems( 0 );
    nitems = fread( &pakeys.SBS, sizeof( unsigned ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LPA_Keys.SBS key from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // read the LPA_db value
    nitems = fread( &pakeys.LPA_db, sizeof( unsigned ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LPA_Keys.LPA_db from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEReader::read( LCI_Keys& cikeys )
  {
    // read the common part of the struct
    LSE_Keys& ekeys = cikeys;
    read( ekeys );

    // read the LCI_script value
    size_t nitems( 0 );
    nitems = fread( &cikeys.LCI_script, sizeof( unsigned ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LCI_Keys LCI_script from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEReader::readKeys( LSE_Keys::KeysType& ktype, LPA_Keys& pakeys, LCI_Keys& cikeys )
  {
    // read the keytype from the file
    size_t nitems(0);
    int itype(0);
    nitems = fread( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
    }
    ktype = static_cast<LSE_Keys::KeysType>( itype );

    // operate on the proper keys object based on the type information
    switch ( ktype ) {
    case LSE_Keys::LPA:
      read( pakeys );
      break;
    case LSE_Keys::LCI:
      read( cikeys );
      break;
    default:
      std::ostringstream ess;
      ess << "LSEReader::read: unknown LSE_Keys typeid " << ktype;
      ess << " from " << m_name;
      throw std::runtime_error( ess.str() );
    }
  }

  void LSEReader::readInfo( LSE_Info::InfoType& infotype, LPA_Info& pinfo, 
			    LCI_ACD_Info& ainfo, LCI_CAL_Info& cinfo, LCI_TKR_Info& tinfo )
  {
    // read in the LSE_Info typeid
    int itype(0);
    size_t nitems(0);
    nitems = fread( &itype, sizeof( int ), 1, m_FILE );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "LSEReader::read: error reading LSE_Info typeid from " << m_name;
      ess << " (" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
    infotype = static_cast<LSE_Info::InfoType>( itype );

    // read the LSE_Info object type id
    unsigned char buf[ 128 * 1024 ];
    size_t len(0);

    // assign to the proper type of object
    switch ( infotype ) {
    case LSE_Info::LPA:
      pinfo.read( m_FILE );
      break;
    case LSE_Info::LCI_ACD:
      read( &buf[0], len );
      ainfo = *reinterpret_cast<LCI_ACD_Info*>( buf );
      break;
    case LSE_Info::LCI_CAL:
      read( &buf[0], len );
      cinfo = *reinterpret_cast<LCI_CAL_Info*>( buf );
      break;
    case LSE_Info::LCI_TKR:
      read( &buf[0], len );
      tinfo = *reinterpret_cast<LCI_TKR_Info*>( buf );
      break;
    default:
      break;
    }
  }

  bool LSEReader::read( LSE_Context& ctx, EBF_Data& ebf, LSE_Info::InfoType& infotype,
			LPA_Info& pinfo, LCI_ACD_Info& ainfo, LCI_CAL_Info& cinfo, LCI_TKR_Info& tinfo,
			LSE_Keys::KeysType& ktype, LPA_Keys& pakeys, LCI_Keys& cikeys )
  {
    // read the context and EBF 
    if ( !read( ctx, ebf ) ) {
      return false;
    }

    // read the LSE_Info object content
    readInfo( infotype, pinfo, ainfo, cinfo, tinfo );

    // read the translated-key objects
    readKeys( ktype, pakeys, cikeys );

    return true;
  }

}
