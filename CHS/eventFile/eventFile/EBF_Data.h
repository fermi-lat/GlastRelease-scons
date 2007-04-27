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

#ifndef EVENTFILE_EBF_DATA_HH
#define EVENTFILE_EBF_DATA_HH

#include <stdio.h>

class EBFevent;

namespace eventFile {

  class EBF_Data {
  public:

    EBF_Data() {};

    const EBFevent* start() const { return reinterpret_cast< const EBFevent* >( &m_data[0] ); };
    const EBFevent* end() const { return reinterpret_cast< const EBFevent* >( &m_data[m_len] ); };
    const unsigned char* data() const { return m_data; }
    unsigned size() const { return m_len; };
    void init( unsigned nbytes, const void* payload );
    void write( FILE* ) const;
    void read( FILE* );

  private:
    unsigned char m_data[128*1024];
    unsigned      m_len;
  };

};

#endif
