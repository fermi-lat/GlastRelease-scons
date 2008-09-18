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
#define LSEHEADER_ALIAS_LEN 64

namespace eventFile {
  
  struct LSEHeader {
    // ctor/dtor
    LSEHeader();
    ~LSEHeader() {};

    // data members
    unsigned m_version;
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

    // version accessor
    unsigned version() const { return m_version & 0x000000FF; }

    // static-variable accessors/mutators
    static unsigned    moot_key()   { return m_moot_key; }
    static const char* moot_alias() { return m_moot_alias; }
    void set_moot_key( unsigned );
    void set_moot_alias( const char* );

    private:

    // static members
    static unsigned m_moot_key;
    static char     m_moot_alias[LSEHEADER_ALIAS_LEN];

    // file-format version specifier
    static const unsigned FormatVersion = 0x09090909;
  };
};

#endif
