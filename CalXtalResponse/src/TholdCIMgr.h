#ifndef TholdCIMgr_H
#define TholdCIMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdCI.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

class TholdCIMgr : public CalibItemMgr {
 public:
  TholdCIMgr() : CalibItemMgr(CalibData::CAL_TholdCI) {};

  /// retrieve threshold calibration constants as measured w/ charge injection
  StatusCode getTholds(const CalXtalId &xtalId,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE,
                       CalibData::ValSig &LAC);

  /// retrieve Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
  StatusCode getULD(const CalXtalId &xtalId,
                    CalibData::ValSig &ULDThold);

  /// retrieve pedestal calibration constants as measured during charge injection threshold testing.
  StatusCode getPed(const CalXtalId &xtalId,
                    CalibData::ValSig &ped);
 private:
  CalibData::CalTholdCI *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::CalTholdCI*) (m_rngBases[XtalIdx(xtalId)]);
  }
  
  bool validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase);
  
  StatusCode fillRangeBases();
  
  bool checkXtalId(const CalXtalId &xtalId) {
    if (!xtalId.validFace())
      throw invalid_argument("TholdCI calib_type requires valid face information."
                             " Programmer error");
    return true;
  }

};

#endif
