#ifndef MPDMgr_H
#define MPDMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalMevPerDac.h"

// EXTLIB

// STD

using namespace CalUtil;
using namespace idents;

using CalibData::ValSig;

class CalCalibSvc;

/// \brief Manage MevPerDAC calibration type
class MPDMgr : public CalibItemMgr {
 public:
  MPDMgr() : 
    CalibItemMgr(CalibData::CAL_MevPerDac) {};

  /// get MeVPerDac ratios for given xtal
  StatusCode getMPD(XtalIdx xtalIdx,
                    ValSig &mpdLrg,
                    ValSig &mpdSm);

 private:
  StatusCode loadIdealVals();

  
  StatusCode genLocalStore();

  bool validateRangeBase(CalibData::CalMevPerDac *mpd);

  ValSig m_idealMPDLrg;
  ValSig m_idealMPDSm;

  bool validateRangeBase(CalibData::RangeBase *rangeBase) {return true;}

};

#endif
