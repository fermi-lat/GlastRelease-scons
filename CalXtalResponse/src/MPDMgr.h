#ifndef MPDMgr_H
#define MPDMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalMevPerDac.h"

// EXTLIB

// STD

using namespace CalDefs;
using namespace idents;

class MPDMgr : public CalibItemMgr {
 public:
  MPDMgr() : CalibItemMgr(CalibData::CAL_MevPerDac) {};

  /// retrieve MeVPerDac ratios for given xtal
  StatusCode getMPD(const CalXtalId &xtalId,
                    CalibData::ValSig &lrg,
                    CalibData::ValSig &sm);

 private:
  CalibData::CalMevPerDac *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::CalMevPerDac*)(m_rngBases[XtalIdx(xtalId)]);
  }
  bool validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase);
  StatusCode fillRangeBases();
  bool checkXtalId(const CalXtalId &xtalId) {return true;}
};

#endif
