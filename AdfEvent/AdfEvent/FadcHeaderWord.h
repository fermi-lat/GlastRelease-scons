#ifndef FADCHEADERWORD_HH
#define FADCHEADERWORD_HH

#include <iostream>
#include "AdfEvent/AncillaryWord.h"
/**
* @class FadcHeaderWord 
* @brief Base class describing a header word from the FADC (Flash Analog to Digital Converter).

* The basic structure is as follows:
* 31             28  27       24  23    12  11      0
* -----------------  -----------  --------  ---------
* ANCILLARY_FADC_ID  fadc_module  (unused)  fadc_fifo
*/

namespace AncillaryData {
  
  class FadcHeaderWord : public AncillaryWord {
  public:
    unsigned int getFadcModule() const {return (m_word >> 24) & 0xf;}
    unsigned int getFadcFifo()   const {return m_word & 0xfff;}
    bool checkHeader()           const {return (getHeader() == ANCILLARY_FADC_ID);}
  };
}//namespace AdfEvent
#endif
