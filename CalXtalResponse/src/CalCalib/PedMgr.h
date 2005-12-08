#ifndef PedMgr_H
#define PedMgr_H

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

class CalCalibSvc;

/** @class PedMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class PedMgr : public CalibItemMgr {
 public:
  PedMgr() : 
    CalibItemMgr(CalibData::CAL_Ped)
    {};

  /// get pedestal vals for given xtal/face/rng
  StatusCode getPed(RngIdx rngIdx,
                    float &avr,
                    float &sig,
                    float &cos);
 private:
  StatusCode loadIdealVals();

  bool validateRangeBase(CalibData::Ped *ped);

  /// ped vals to use when calib db is down
  CalArray<RngNum,float> m_idealPeds;   
  /// ped sigma vals to use when calib db is down
  CalArray<RngNum,float> m_idealPedSig; 
  /// correlated ped cosine vals to use when calib db is down
  CalArray<RngNum,float> m_idealCos;    

  StatusCode genLocalStore();
};

#endif
