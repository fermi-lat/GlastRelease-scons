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

using CalibData::ValSig;

class MPDMgr : public CalibItemMgr {
 public:
  MPDMgr(const IdealCalCalib &idealCalib) : 
    CalibItemMgr(CalibData::CAL_MevPerDac,
                 idealCalib) {};

  /// retrieve MeVPerDac ratios for given xtal
  StatusCode getMPD(const CalXtalId &xtalId,
                    ValSig &mpdLrg,
                    ValSig &mpdSm);

 private:
  CalibData::CalMevPerDac *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::CalMevPerDac*)(m_rngBases[XtalIdx(xtalId)]);
  }

  bool validateRangeBase(const CalXtalId&, CalibData::RangeBase*)
    {return true;}

  StatusCode fillRangeBases();

  bool checkXtalId(const CalXtalId&) {return true;}

  StatusCode loadIdealVals();

  ValSig m_idealMPDLrg;
  ValSig m_idealMPDSm;
};

#endif
