#ifndef ANCILLARYDATATAIL_HH
#define ANCILLARYDATATAIL_HH

#include <iostream>

/**
* @class AncillaryDataTail
*/

const unsigned int ANCILLARY_TAIL_ID      = 14;

namespace AncillaryData {
  
  class AncillaryDataTail {
  public:
    AncillaryDataTail(){;} 
    ~AncillaryDataTail()                    {;} 
    void setData(unsigned int word[3]) {m_word = word;}
  private:
    unsigned int *m_word;
  };
  
}
#endif

