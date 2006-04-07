/**
 * @class eventFile::LSEReader
 *
 * @brief Class for reading a file containing per-event context, meta-info, and EBF data
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSEREADER_H
#define LSEREADER_H

#include <stdio.h>

#include <string>

#include "eventFile/LSE_Info.h"

namespace eventFile {

  class LSE_Context;
  class EBF_Data;
  
  class LSEReader {
  public:
    LSEReader( const std::string& filename );
    ~LSEReader();

    bool read( LSE_Context&, EBF_Data&, LSE_Info::InfoType&, LPA_Info&, LCI_ACD_Info&, LCI_CAL_Info&, LCI_TKR_Info& );
    void close();
    unsigned runid() const { return m_runid; };

    unsigned long long evtcnt() const { return m_evtcnt; };
    unsigned long long begGEM() const { return m_GEMseq_beg; };
    unsigned long long endGEM() const { return m_GEMseq_end; };

  private:
    std::string m_name;
    unsigned m_runid;
    unsigned long long m_evtcnt;
    unsigned long long m_GEMseq_beg;
    unsigned long long m_GEMseq_end;
    FILE* m_FILE;

    bool read( LSE_Context&, EBF_Data& );
    void read( int&, unsigned char*, size_t& );
    void readHeader();
  };
  
};

#endif
