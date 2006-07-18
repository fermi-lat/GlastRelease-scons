#ifndef SCALERHEADERWORD_HH
#define SCALERHEADERWORD_HH

#include <iostream>
#include "AdfEvent/AncillaryWord.h"
/**
* @class FadcHeaderWord 
* @brief Base class describing a header word from the scaler.

* The basic structure is as follows:
* 31               28  27     8  7         0
* -------------------  --------  -----------
* ANCILLARY_SCALER_ID  (unused)  scaler_fifo
*/

namespace AncillaryData {

	class ScalerHeaderWord : public AncillaryWord {
	public:
	  unsigned int getScalerCounters() {return m_word & 0xff;}
	  bool checkHeader()           const {return (getHeader() == ANCILLARY_SCALER_ID);}
	};
	
}//namespace AdfEvent
#endif
