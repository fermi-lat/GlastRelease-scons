/** @class eventFile::EBF_Data
 *
 * @brief Class encapsulating the EBF representation of an event.
 *
 * This class encapsulates the EBF representation of an event, as unpacked / uncompressed
 * from the downlink representation by the DFI package.  This representaion is suitable
 * for parsing by the LDF library.  The start() and end() methods provide pointers delineating
 * the range of memory over which an EBFeventIterator subclass should operate.
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef EBF_DATA_HH
#define EBF_DATA_HH

class EBFevent;

namespace eventFile {

  class EBF_Data {
  public:

    EBF_Data() {};

    const EBFevent* start() const { return reinterpret_cast< const EBFevent* >( &data[0] ); };
    const EBFevent* end() const { return reinterpret_cast< const EBFevent* >( &data[len] ); };
    unsigned size() const { return len; };

  private:
    unsigned char data[128*1024];
    unsigned      len;
    void init( unsigned nbytes, const void* payload ) {
      *( reinterpret_cast< int* >( &data[0] ) ) = 0x104f0010  ;
      *( reinterpret_cast< int* >( &data[4] ) ) = nbytes +  8 ;
      memcpy( &data[8], payload, nbytes );
      len = nbytes + 8;
    };

  };

};

#endif
