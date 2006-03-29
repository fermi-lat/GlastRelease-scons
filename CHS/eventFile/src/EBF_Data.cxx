#include <errno.h>
#include <stdexcept>
#include <sstream>

#include "eventFile/EBF_Data.h"

namespace eventFile {

  void EBF_Data::write( FILE* fp ) const
  {
    size_t nitems(0);

    // write out the size of the EBF blob
    nitems = fwrite( &len, sizeof( len ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "EBF_Data::write: error writing length ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // write out the ebf blob itself
    nitems = fwrite( data, len, 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "EBF_Data::write: error writing data ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void EBF_Data::read( FILE* fp )
  {
    size_t nitems(0);

    // read in the length of the EBF blob
    nitems = fread( &len, sizeof( len ), 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "EBF_Data::read: error reading length ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }

    // read in the EBF blob itself
    nitems = fread( data, len, 1, fp );
    if ( nitems != 1 ) {
      std::ostringstream ess;
      ess << "EBF_Data::read: error reading length ";
      ess << "(" << errno << "=" << strerror( errno ) << ")";
      throw std::runtime_error( ess.str() );
    }
  }

  void EBF_Data::init( unsigned nbytes, const void* payload )
  {
    *( reinterpret_cast< int* >( &data[0] ) ) = 0x104f0010  ;
    *( reinterpret_cast< int* >( &data[4] ) ) = nbytes +  8 ;
    memcpy( &data[8], payload, nbytes );
    len = nbytes + 8;
  }

}
