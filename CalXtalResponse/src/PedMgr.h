#ifndef PedMgr_H
#define PedMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/Ped.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

class PedMgr : public CalibItemMgr {
 public:
  PedMgr() : CalibItemMgr(CalibData::CAL_Ped) {};

  /// retrieve pedestal values for given xtal/face/rng
  StatusCode getPed(const CalXtalId &xtalId,
                    float &avr,
                    float &sig,
                    float &cos);
 private:
  CalibData::Ped *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::Ped*)(m_rngBases[RngIdx(xtalId)]);
  }
  
  bool validateRangeBase(const CalXtalId&, CalibData::RangeBase*)
    {return true;}
  
  StatusCode fillRangeBases();
  
  bool checkXtalId(const CalXtalId &xtalId) {
    if (!xtalId.validRange() || !xtalId.validFace())
      throw invalid_argument("Pedestal calib_type requires valid range & face info in CalXtalId");
    return true;
  }
};

#endif
