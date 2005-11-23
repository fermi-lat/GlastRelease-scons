#ifndef TholdCIMgr_H
#define TholdCIMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdCI.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

using CalibData::ValSig;

class CalCalibSvc;

/** @class TholdCIMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal charge-injection measured threshold data.
*/

class TholdCIMgr : public CalibItemMgr {
 public:
  TholdCIMgr() 
    : CalibItemMgr(CalibData::CAL_TholdCI),
    m_idealULD(RngNum::N_VALS),
    m_idealPed(RngNum::N_VALS) {};

  /// get threshold calibration constants as measured w/ charge injection
  StatusCode getTholds(CalXtalId xtalId,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE,
                       CalibData::ValSig &LAC);

  /// get Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
  StatusCode getULD(CalXtalId xtalId,
                    CalibData::ValSig &ULDThold);

  /// get pedestal calibration constants as measured during charge injection threshold testing.
  StatusCode getPed(CalXtalId xtalId,
                    CalibData::ValSig &ped);
 private:
  bool checkXtalId(CalXtalId xtalId);

  StatusCode loadIdealVals();

  LATWideIndex genIdx(CalXtalId xtalId) {return FaceIdx(xtalId);}

  ValSig m_idealFLE;
  ValSig m_idealFHE;
  ValSig m_idealLAC;
  CalVec<RngNum, ValSig> m_idealULD;
  CalVec<RngNum, ValSig> m_idealPed;

  bool validateRangeBase(CalibData::RangeBase *rangeBase);
  
};

#endif
