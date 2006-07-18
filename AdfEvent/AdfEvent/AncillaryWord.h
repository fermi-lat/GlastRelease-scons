#ifndef ANCILLARYWORD_HH
#define ANCILLARYWORD_HH
#include "AdfEvent/TaggerParameters.h"

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
  const unsigned int ANCILLARY_WORD_LENGTH  = 4;

  const unsigned int ANCILLARY_QDC_ID       = 1;
  const unsigned int ANCILLARY_QDC_HID      = 3;

  const unsigned int ANCILLARY_FADC_ID      = 2;
  const unsigned int ANCILLARY_FADC_HID     = 4;

  const unsigned int ANCILLARY_SCALER_ID    = 10;

  const unsigned int ANCILLARY_END_RUN      = 14;


  class AncillaryWord {
  public:
    AncillaryWord() {;}
    AncillaryWord(unsigned int word) {m_word = word;}
    //    AncillaryWord(const AncillaryWord &aword) {m_word=aword.getRawWord();}
    void setData(unsigned int word) {m_word = word;}
    ~AncillaryWord()                 {;} 
    unsigned int  getRawWord()  const {return m_word;}
    unsigned int  getHeader()   const {return (m_word >> 28) & 0xf;}
    unsigned int  getRawData()  const {return m_word & 0xfffffff;}
    
  protected:
    unsigned int m_word;
    
  };
  
}//namespace AdfEvent
#endif

