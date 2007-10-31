#ifndef LSFREADER_H
#define LSFREADER_H

#include "eventFile/LSEReader.h"

namespace eventFile {

  class LSE_Context;
  class EBF_Data;

};

namespace lsfData {
  class LsfCcsds;
  class MetaEvent;
  class LciConfiguration;

  class LSFReader : public eventFile::LSEReader {
  public:
    LSFReader( const std::string& filename ) : eventFile::LSEReader( filename ) {};
    ~LSFReader() {};

    bool read( LsfCcsds&, MetaEvent&, eventFile::EBF_Data& );

    void transferCcsds( const eventFile::LSE_Context&, LsfCcsds& );
    void transferContext( const eventFile::LSE_Context&, MetaEvent& );
    void transferTime( const eventFile::LSE_Context&, const eventFile::LSE_Info&,     MetaEvent& );
    void transferLciCfg( const eventFile::LCI_Info&, LciConfiguration& );
    void transferInfo( const eventFile::LSE_Context&, const eventFile::LPA_Info&,     MetaEvent& );
    void transferInfo( const eventFile::LSE_Context&, const eventFile::LCI_ACD_Info&, MetaEvent& );
    void transferInfo( const eventFile::LSE_Context&, const eventFile::LCI_CAL_Info&, MetaEvent& );
    void transferInfo( const eventFile::LSE_Context&, const eventFile::LCI_TKR_Info&, MetaEvent& );
    void transferKeys( const eventFile::LPA_Keys&, MetaEvent& );
    void transferKeys( const eventFile::LCI_Keys&, MetaEvent& );
  };
};

#endif
