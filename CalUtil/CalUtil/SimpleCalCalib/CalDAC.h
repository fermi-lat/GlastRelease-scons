#ifndef CalDAC_h
#define CalDAC_h

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

  /** \brief Represents GLAST Cal DAC threhsold settings (one integer per xtal face)

  contains read & write methods to various file formats

  @author Zachary Fewtrell
  */
  class CalDAC {
  public:
    CalDAC();

    /// write DACestals to columnar TXTfile
    void writeTXT(const std::string &filename) const;

    /// read DACestals from columnar TXTfile
    void readTXT(const std::string &filename);

    int getDAC(const CalUtil::FaceIdx faceIdx) const {
      return m_DACs[faceIdx];
    }

    void setDAC(const CalUtil::FaceIdx faceIdx,
                const int val) {
      m_DACs[faceIdx] = val;
    }

    /// unset channels are initailized to this value.
    static const int INVALID_DAC;

  private:
    /// DAC data
    CalUtil::CalVec<CalUtil::FaceIdx, int> m_DACs;
  };

}; // namespace CalUtil
#endif
