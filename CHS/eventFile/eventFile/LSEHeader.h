/**
 * @class eventFile::LSEHeader
 *
 * @brief Class representing the header of an LSE event file
 *
 * @author Bryson Lee <blee@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSEHEADER_H
#define LSEHEADER_H

#include <stdio.h>

#define LSEHEADER_MAX_APIDS 4

namespace eventFile {
  
  struct LSEHeader {
    // ctor/dtor
    LSEHeader();
    ~LSEHeader() {};

    // data members
    unsigned m_runid;
    unsigned m_secs_beg;
    unsigned m_secs_end;
    unsigned long long m_evtcnt;
    unsigned long long m_GEMseq_beg;
    unsigned long long m_GEMseq_end;
    unsigned m_src_apids[LSEHEADER_MAX_APIDS];
    unsigned m_src_seqerr[LSEHEADER_MAX_APIDS];
    unsigned m_dfi_apids[LSEHEADER_MAX_APIDS];
    unsigned m_dfi_dfierr[LSEHEADER_MAX_APIDS];

    // read-write routines
    void read( FILE* );
    void write( FILE* );
  };
};

#endif
