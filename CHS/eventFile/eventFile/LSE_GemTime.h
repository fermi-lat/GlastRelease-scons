/**
 * @class eventFile::LSE_GemTime
 *
 * @brief Class encapsulating timestamp information from the GEM
 *
 * This class exists separately in the public interface simply to allow
 * reuse between the LSE_Context and LSE_Info objects.
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSE_GEMTIME_HH
#define LSE_GEMTIME_HH

namespace eventFile {
  
  struct LSE_GemTime {
    LSE_GemTime() {};
    void dump( const char* pre, const char* post ) const;
    unsigned tics;
    unsigned hacks;
  };

};

#endif
