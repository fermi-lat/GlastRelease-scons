/**
 * @class eventFile::LPA_File
 *
 * @brief Class encapsulating a file containing merged LPA data
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LPA_FILE_HH
#define LPA_FILE_HH

#include <stdio.h>

#include <string>

namespace eventFile {

  class LSE_Context;
  class LPA_Info;
  class EBF_Data;
  
  class LPA_File {
  public:
    enum Mode { Read = 0, Write };
    LPA_File( const std::string& filename, Mode openmode, unsigned runid = 0 );
    ~LPA_File();

    void write( const LSE_Context&, const LPA_Info&, const EBF_Data& );
    bool read( LSE_Context&, LPA_Info&, EBF_Data& );
    void close();
    unsigned runid() const { return m_runid; };

    void evtcnt( unsigned long long evt ) { m_evtcnt = evt; };
    void begGEM( unsigned long long gem ) { m_GEMseq_beg = gem; };
    void endGEM( unsigned long long gem ) { m_GEMseq_end = gem; };
    unsigned long long evtcnt() const { return m_evtcnt; };
    unsigned long long begGEM() const { return m_GEMseq_beg; };
    unsigned long long endGEM() const { return m_GEMseq_end; };

  private:
    std::string m_name;
    std::string m_mode;
    unsigned m_runid;
    unsigned long long m_evtcnt;
    unsigned long long m_GEMseq_beg;
    unsigned long long m_GEMseq_end;
    FILE* m_FILE;

  };
  
};

#endif
