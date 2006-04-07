#ifndef LSFREADER_H
#define LSFREADER_H

#include "eventFile/LSEReader.h"

namespace lsfData {
  class LsfCcsds;
  class MetaEvent;
  class LciConfiguration;
};

namespace eventFile {

  class LSE_Context;
  class EBF_Data;

  class LSFReader : public LSEReader {
  public:
    LSFReader( const std::string& filename ) : LSEReader( filename ) {};
    ~LSFReader() {};

    bool read( lsfData::LsfCcsds&, lsfData::MetaEvent&, EBF_Data& );

    void transferCcsds( const LSE_Context&, lsfData::LsfCcsds& );
    void transferContext( const LSE_Context&, lsfData::MetaEvent& );
    void transferTime( const LSE_Context&, const LSE_Info&,     lsfData::MetaEvent& );
    void transferLciCfg( const LCI_Info&, lsfData::LciConfiguration& );
    void transferInfo( const LSE_Context&, const LPA_Info&,     lsfData::MetaEvent& );
    void transferInfo( const LSE_Context&, const LCI_ACD_Info&, lsfData::MetaEvent& );
    void transferInfo( const LSE_Context&, const LCI_CAL_Info&, lsfData::MetaEvent& );
    void transferInfo( const LSE_Context&, const LCI_TKR_Info&, lsfData::MetaEvent& );
  };
};

#endif
