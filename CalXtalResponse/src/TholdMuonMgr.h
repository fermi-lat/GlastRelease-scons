#ifndef TholdMuonMgr_H
#define TholdMuonMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdMuon.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

class TholdMuonMgr : public CalibItemMgr {
 public:
  TholdMuonMgr() : CalibItemMgr(CalibData::CAL_TholdMuon) {};

  /// retrieve threshold calibration constants as measured w/ muon calibration
  StatusCode getTholds(const CalXtalId &xtalId,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE);

  /// retrieve pedestal calibration constants as measured during muon calibration threshold testing.
  StatusCode getPed(const CalXtalId &xtalId,
                    CalibData::ValSig &ped);
 private:
  CalibData::CalTholdMuon *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::CalTholdMuon*)(m_rngBases[XtalIdx(xtalId)]);
  }
  
  bool validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase);
  
  StatusCode fillRangeBases();
  
  bool checkXtalId(const CalXtalId &xtalId) {
    if (!xtalId.validFace())
      throw invalid_argument("TholdMuon calib_type requires valid face information."
                             " Programmer error");
    return true;
  }

};

#endif
