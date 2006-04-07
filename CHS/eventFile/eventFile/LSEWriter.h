/**
 * @class eventFile::LSEWriter
 *
 * @brief Class for writing per-event context, meta-info, and EBF data
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSEWRITER_H
#define LSEWRITER_H

#include <stdio.h>

#include <string>
#include <utility>

#include "eventFile/LSEHeader.h"

namespace eventFile {

  class LSE_Context;
  class EBF_Data;
  class LPA_Info;
  class LCI_ACD_Info;
  class LCI_CAL_Info;
  class LCI_TKR_Info;
  
  class LSEWriter {
  public:
    LSEWriter( const std::string& filename, unsigned runid = 0 );
    ~LSEWriter();

    void write( const LSE_Context&, const EBF_Data&, const LPA_Info& );
    void write( const LSE_Context&, const EBF_Data&, const LCI_ACD_Info& );
    void write( const LSE_Context&, const EBF_Data&, const LCI_CAL_Info& );
    void write( const LSE_Context&, const EBF_Data&, const LCI_TKR_Info& );

    void close();

    // header mutators
    void runid( unsigned rid ) { m_hdr.m_runid = rid; };
    void begSec( unsigned sec ) { m_hdr.m_secs_beg = sec; };
    void endSec( unsigned sec ) { m_hdr.m_secs_end = sec; };
    void evtcnt( unsigned long long evt ) { m_hdr.m_evtcnt = evt; };
    void begGEM( unsigned long long gem ) { m_hdr.m_GEMseq_beg = gem; };
    void endGEM( unsigned long long gem ) { m_hdr.m_GEMseq_end = gem; };
    void seqErr( unsigned apid, unsigned seqerr, int islot )
      {
	if ( islot >= LSEHEADER_MAX_APIDS ) return;
	m_hdr.m_src_apids[islot]  = apid;
	m_hdr.m_src_seqerr[islot] = seqerr;
      };

    // header accessors
    unsigned runid() const { return m_hdr.m_runid; };
    unsigned begSec() const { return m_hdr.m_secs_beg; };
    unsigned endSec() const { return m_hdr.m_secs_end; };
    unsigned long long evtcnt() const { return m_hdr.m_evtcnt; };
    unsigned long long begGEM() const { return m_hdr.m_GEMseq_beg; };
    unsigned long long endGEM() const { return m_hdr.m_GEMseq_end; };
    std::pair<unsigned, unsigned> seqErr( int islot ) const
      {
	if ( islot < LSEHEADER_MAX_APIDS ) {
	  return std::pair<unsigned, unsigned>( m_hdr.m_src_apids[islot], m_hdr.m_src_seqerr[islot] );
	} else {
	  return std::pair<unsigned, unsigned>( 0, 0 );
	}
      };

  private:
    std::string m_name;
    LSEHeader m_hdr;
    FILE* m_FILE;

    void write( const LSE_Context&, const EBF_Data& );
    void write( int, const void*, size_t );
    void writeHeader();
  };
  
};

#endif
