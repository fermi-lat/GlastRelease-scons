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

  private:
    std::string m_name;
    std::string m_mode;
    unsigned m_runid;
    FILE* m_FILE;

  };
  
};

#endif
