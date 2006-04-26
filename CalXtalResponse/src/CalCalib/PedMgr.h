#ifndef PedMgr_H
#define PedMgr_H
// $Header$
// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/Ped.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB
// STD

using namespace CalUtil;
using namespace idents;
using namespace CalibData;

class CalCalibSvc;

/** @class PedMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class PedMgr : public CalibItemMgr {
 public:
  PedMgr(CalCalibShared &ccsShared) : 
    CalibItemMgr(CAL_Ped, ccsShared)
    {};

  /// get pedestal vals for given xtal/face/rng
  const Ped *getPed(RngIdx rngIdx);
 private:
  StatusCode loadIdealVals();

  bool validateRangeBase(Ped *ped);
  
  CalArray<RngNum, Ped> m_idealPeds;

  StatusCode genLocalStore();
};

#endif
