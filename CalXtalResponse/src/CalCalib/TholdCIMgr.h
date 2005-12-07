#ifndef TholdCIMgr_H
#define TholdCIMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/CalTholdCI.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB
// STD

using namespace CalUtil;
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
    : CalibItemMgr(CalibData::CAL_TholdCI)
    {};

  /// get threshold calibration constants as measured w/ charge injection
  StatusCode getTholds(FaceIdx faceIdx,
                       CalibData::ValSig &FLE,
                       CalibData::ValSig &FHE,
                       CalibData::ValSig &LAC);

  /// get Upper Level Discriminator threshold as measured w/ charnge
  /// injection for given xtal/face/rng

  StatusCode getULD(RngIdx rngIdx,
                    CalibData::ValSig &ULDThold);

  /// get pedestal calibration constants as measured during charge
  /// injection threshold testing.

  StatusCode getPed(RngIdx rngIdx,
                    CalibData::ValSig &ped);
 private:
  StatusCode loadIdealVals();

  
  StatusCode genLocalStore();

  ValSig m_idealFLE;
  ValSig m_idealFHE;
  ValSig m_idealLAC;
<<<<<<< TholdCIMgr.h
  CalArray<RngNum, ValSig> m_idealULD;
  CalArray<RngNum, ValSig> m_idealPed;

  
  /// Validate TDS data entry (for empty ptrs & fun stuff like that)
  bool validateRangeBase(CalibData::CalTholdCI *tholdCI);
=======
  CalVec<RngNum, ValSig> m_idealULD;
  CalVec<RngNum, ValSig> m_idealPed;

  bool validateRangeBase(CalibData::RangeBase *rangeBase);
>>>>>>> 1.4
  
};

#endif
