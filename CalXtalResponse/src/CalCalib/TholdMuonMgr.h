#ifndef TholdMuonMgr_H
#define TholdMuonMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdMuon.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB
// STD

using namespace CalUtil;
using namespace idents;

using CalibData::ValSig;

class CalCalibSvc;

/** @class TholdMuonMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal muon-measured threshold calibration data.
*/

class TholdMuonMgr : public CalibItemMgr {
 public:
  TholdMuonMgr() : 
    CalibItemMgr(CalibData::CAL_TholdMuon)
    {};

  /// get threshold calibration constants as measured w/ muon calibration
  StatusCode getTholds(FaceIdx faceIdx,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE);

  /// get pedestal calibration constants as measured during muon calibration threshold testing.
  StatusCode getPed(RngIdx rngIdx,
                    CalibData::ValSig &ped);
 private:
  StatusCode loadIdealVals();

  
  StatusCode genLocalStore();

  ValSig m_idealFLE;
  ValSig m_idealFHE;
  CalArray<RngNum, ValSig> m_idealPed;

  /// Validate TDS data entry (for empty ptrs & fun stuff like that)
  bool validateRangeBase(CalibData::CalTholdMuon *tholdMuon); 

  bool validateRangeBase(CalibData::RangeBase *rangeBase);
};

#endif
