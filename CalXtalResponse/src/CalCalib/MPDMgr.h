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

class CalCalibSvc;

class MPDMgr : public CalibItemMgr {
 public:
  MPDMgr() : 
    CalibItemMgr(CalibData::CAL_MevPerDac) {};

  /// get MeVPerDac ratios for given xtal
  StatusCode getMPD(const CalXtalId &xtalId,
                    ValSig &mpdLrg,
                    ValSig &mpdSm);

 private:
  bool validateRangeBase(CalibData::RangeBase*);

  StatusCode fillRangeBases();

  bool checkXtalId(const CalXtalId&) {return true;}

  StatusCode loadIdealVals();

  LATWideIndex genIdx(const CalXtalId &xtalId) {return XtalIdx(xtalId);}

  ValSig m_idealMPDLrg;
  ValSig m_idealMPDSm;
};

#endif
