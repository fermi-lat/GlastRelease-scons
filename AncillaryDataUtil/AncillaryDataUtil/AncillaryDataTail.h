#ifndef ANCILLARYDATATAIL_HH
#define ANCILLARYDATATAIL_HH

#include <iostream>

/**
* @class AncillaryDataTail
*/

namespace AncillaryData {
  const unsigned int ANCILLARY_TAIL_VERSION = 100;
  const unsigned int ANCILLARY_TAIL_ID      = 14;
  const unsigned int TAIL_LENGTH            =5; 
 
  class AncillaryDataTail {
  public:
    AncillaryDataTail(){;} 
    ~AncillaryDataTail()                    {;} 
    void setData(unsigned int word[TAIL_LENGTH]) {m_word = word;}
    unsigned int totalNumEvents(){return m_word[2] & 0xfffffff;}
  private:
    unsigned int *m_word;
    
  };
  
}
#endif

