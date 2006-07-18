#ifndef ANCILLARYDATAHEADER_HH
#define ANCILLARYDATAHEADER_HH

#include <iostream>

/**
* @class AncillaryDataHeader
*/

namespace AncillaryData {
  
  const unsigned int ANCILLARY_HEADER_VERSION = 200;
  const unsigned int ANCILLARY_FILE_HEADER  = 0xf1040;
  const unsigned int ANCILLARY_HEADER_ID    = 15;
  const unsigned int HEADER_LENGTH          = 6; 
  
  class AncillaryDataHeader {
  public:
    AncillaryDataHeader(){;} 
    ~AncillaryDataHeader()                    {;} 
    void setData(unsigned *word) 
      {
	m_version       = unsigned((word[0] >> 20) & 0xfff);
	m_length        = int(word[1]);
	m_runNumber     = int(word[2]);
	m_commentLength = int(word[5] & 0xfffffff);
      }
    unsigned getVersion()        const {return m_version;}
    unsigned getHeaderLength()   const {return m_length;}
    unsigned getRunNumber()      const {return m_runNumber;}
    unsigned getCommentLength()  const {return m_commentLength;}
    
    bool checkVersion()              const {return (getVersion() ==  ANCILLARY_HEADER_VERSION );}
    
  private:
    unsigned m_version;
    unsigned m_length;
    unsigned m_runNumber;
    unsigned m_commentLength;
  };
  
}//namespace AdfEvent
#endif

