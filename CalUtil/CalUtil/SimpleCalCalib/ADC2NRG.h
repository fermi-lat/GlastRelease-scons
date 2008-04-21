#ifndef ADC2NRG_h
#define ADC2NRG_h

// $Header$

/** @file
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// STD INCLUDES

namespace CalUtil {

  /** \brief Represents GLAST Cal adc2mev conversion factor for each ADC range

  contains read & write methods to various file formats

  @author Zachary Fewtrell
  */
  class ADC2NRG {
  public:
    ADC2NRG();

    /// write DACestals to columnar TXTfile
    void writeTXT(const std::string &filename) const;

    /// read DACestals from columnar TXTfile
    void readTXT(const std::string &filename);

    float getADC2NRG(const CalUtil::RngIdx rngIdx) const {
      return m_adc2nrg[rngIdx];
    }

    void setADC2NRG(const CalUtil::RngIdx rngIdx,
                const int val) {
      m_adc2nrg[rngIdx] = val;
    }

    /// unset channels are initailized to this value.
    static const float INVALID_ADC2NRG;

  private:
    /// DAC data
    CalUtil::CalVec<CalUtil::RngIdx, float> m_adc2nrg;
  };

}; // namespace CalUtil
#endif
