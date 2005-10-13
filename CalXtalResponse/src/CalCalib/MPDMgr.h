#ifndef MPDMgr_H
#define MPDMgr_H

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
  StatusCode getMPD(CalXtalId xtalId,
                    ValSig &mpdLrg,
                    ValSig &mpdSm);

 private:
  bool checkXtalId(CalXtalId);

  StatusCode loadIdealVals();

  LATWideIndex genIdx(CalXtalId xtalId) {return XtalIdx(xtalId);}

  ValSig m_idealMPDLrg;
  ValSig m_idealMPDSm;
};

#endif
