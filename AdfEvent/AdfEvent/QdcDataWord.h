#ifndef QDCDATAWORD_HH
#define QDCDATAWORD_HH

#include <iostream>
#include "AdfEvent/AncillaryWord.h"
/**
* @class QdcDataWord 
* @brief Base class describing a data word from the QDC (Q to Digital Converter).

* The basic structure is as follows:
* 31            28  27    20  19       12  11      0
* ----------------  --------  -----------  ---------
* ANCILLARY_QDC_ID  (unused)  qdc_channel  qdc_value
*/

namespace AncillaryData {
  class QdcDataWord : public AncillaryWord {
  public:
    unsigned int getQdcModule()  const {return 0;}
    unsigned int getQdcChannel() const {return (m_word >> 12) & 0xff;}
    unsigned int getQdcValue()   const {return m_word & 0xfff;}
    bool checkHeader()           const {return (getHeader() == ANCILLARY_QDC_ID);}
  };
}//namespace AdfEvent
#endif
  
