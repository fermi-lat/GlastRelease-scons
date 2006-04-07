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

    void evtcnt( unsigned long long evt ) { m_evtcnt = evt; };
    void begGEM( unsigned long long gem ) { m_GEMseq_beg = gem; };
    void endGEM( unsigned long long gem ) { m_GEMseq_end = gem; };
    unsigned long long evtcnt() const { return m_evtcnt; };
    unsigned long long begGEM() const { return m_GEMseq_beg; };
    unsigned long long endGEM() const { return m_GEMseq_end; };
    unsigned runid() const { return m_runid; };


  private:
    std::string m_name;
    unsigned m_runid;
    unsigned long long m_evtcnt;
    unsigned long long m_GEMseq_beg;
    unsigned long long m_GEMseq_end;
    FILE* m_FILE;

    void write( const LSE_Context&, const EBF_Data& );
    void write( int, const void*, size_t );
    void writeHeader();
  };
  
};

#endif
