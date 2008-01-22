#ifndef TholdCI_h
#define TholdCI_h

/** @file 
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// GLAST INCLUDES

// EXTLIB INCLUDES

// STD INCLUDES
#include <string>

namespace CalUtil {

  /** represent full Cal_TholdCI calibration data set
   */
  class TholdCI {
  public:
    TholdCI();

    /// read thresholds from txt file
    void readTXT(const std::string &filename);

    float getFLEThresh(const FaceIdx faceIdx) const {
      return m_fleThresh[faceIdx];
    }

    float getFHEThresh(const FaceIdx faceIdx) const {
      return m_fheThresh[faceIdx];
    }

    float getLACThresh(const FaceIdx faceIdx) const {
      return m_lacThresh[faceIdx];
    }

    float getULDThresh(const RngIdx rngIdx) const {
      return m_uldThresh[rngIdx];
    }

    /// used as a placeholder for unknown calibrations
    static const short INVALID_THOLD;

  private:
    CalVec<FaceIdx, float> m_fleThresh;
    CalVec<FaceIdx, float> m_fheThresh;
    CalVec<FaceIdx, float> m_lacThresh;
    CalVec<RngIdx, float> m_uldThresh;
    CalVec<RngIdx, float> m_ped;

  };
} 

#endif // TholdCI_h
