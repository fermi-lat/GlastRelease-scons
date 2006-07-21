#ifndef FADCDATAWORD_HH
#define FADCDATAWORD_HH

#include <iostream>
#include "AdfEvent/AncillaryWord.h"
#include "AdfEvent/TaggerParameters.h"
/**
* @class FadcDataWord 
* @brief Base class describing a data word from the FADC (Flash Analog to Digital Converter).

* The basic structure is as follows:
* 31             28  27       24  23        12  11       0
* -----------------  -----------  ------------  ----------
* ANCILLARY_FADC_ID  fadc_module  fadc_channel  fadc_value
*/

namespace AncillaryData {
  class FadcDataWord : public AncillaryWord {
  public:
    unsigned int getFadcModule()    const {return (m_word >> 24) & 0xf;}
    unsigned int getFadcChannel()   const {return (m_word >> 12) & 0xfff;}
    unsigned int getTaggerLayer()   const {return ((getFadcChannel()< N_CHANNELS_PER_LAYER) ? 0 :1 );}
    unsigned int getTaggerChannel() const {return getFadcChannel()-getTaggerLayer()*N_CHANNELS_PER_LAYER;}
    unsigned int getFadcValue()   const {return m_word & 0xfff;}
    bool checkHeader()            const {return (getHeader() == ANCILLARY_FADC_ID);}
  };
}//namespace AdfEvent
#endif
