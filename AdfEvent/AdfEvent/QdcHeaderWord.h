#ifndef QDCHEADERWORD_HH
#define QDCHEADERWORD_HH

#include <iostream>
#include "AdfEvent/AncillaryWord.h"
/**
* @class QdcHeaderWord 
* @brief Base class describing a header word from the QDC (Q to Digital Converter).

* The basic structure is as follows:
* 31            28  27     5  4      0
* ----------------  --------  --------
* ANCILLARY_QDC_ID  (unused)  qdc_fifo
*/

namespace AncillaryData {
  class QdcHeaderWord : public AncillaryWord {
  public:
    unsigned int getQcdFifo() const {return m_word & 0x1f;}
    bool checkHeader()        const {return (getHeader() == ANCILLARY_QDC_ID);}
  };
  
}//namespace AdfEvent
#endif
