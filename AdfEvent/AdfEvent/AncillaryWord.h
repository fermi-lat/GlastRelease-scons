#ifndef ANCILLARYWORD_HH
#define ANCILLARYWORD_HH

#include <iostream>

/**
* @class AncillaryWord 
* @brief Base class describing a data word from the ancillary detectors.

* Based on ancillary data version 2.00. The basic structure is as follows:
* 31  28  27   0
* ------  ------
* header   data
* Essentially it consists of a 4 bits header and a data content.
* The structure of this data content is specified, case by case, in several other
* classes inheriting by AncillaryWord.
*/

namespace AncillaryData {
  const unsigned ANCILLARY_WORD_LENGTH  = 4;
  const unsigned ANCILLARY_QDC_ID       = 1;
  const unsigned ANCILLARY_FADC_ID      = 2;
  const unsigned ANCILLARY_END_RUN      = 14;
  const unsigned ANCILLARY_SCALER_ID    = 10;

  class AncillaryWord {
  public:
    AncillaryWord() {;}
    AncillaryWord(unsigned word) {m_word = word;}
    //    AncillaryWord(const AncillaryWord &aword) {m_word=aword.getRawWord();}
    void setData(unsigned word) {m_word = word;}
    ~AncillaryWord()                 {;} 
    unsigned  getRawWord()  const {return m_word;}
    unsigned  getHeader()   const {return (m_word >> 28) & 0xf;}
    unsigned  getRawData()  const {return m_word & 0xfffffff;}
    
  protected:
    unsigned m_word;
    
  };
  
}//namespace AdfEvent
#endif

