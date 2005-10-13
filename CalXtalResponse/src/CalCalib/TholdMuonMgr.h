#ifndef TholdMuonMgr_H
#define TholdMuonMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdMuon.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
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
    CalibItemMgr(CalibData::CAL_TholdMuon),
    m_idealPed(RngNum::N_VALS)
    {};

  /// get threshold calibration constants as measured w/ muon calibration
  StatusCode getTholds(CalXtalId xtalId,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE);

  /// get pedestal calibration constants as measured during muon calibration threshold testing.
  StatusCode getPed(CalXtalId xtalId,
                    CalibData::ValSig &ped);
 private:
  bool checkXtalId(CalXtalId xtalId);

  StatusCode loadIdealVals();

  LATWideIndex genIdx(CalXtalId xtalId) {return FaceIdx(xtalId);}

  ValSig m_idealFLE;
  ValSig m_idealFHE;
  CalVec<RngNum, ValSig> m_idealPed;

};

#endif
