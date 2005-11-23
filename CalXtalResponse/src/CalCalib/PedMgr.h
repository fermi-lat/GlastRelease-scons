#ifndef PedMgr_H
#define PedMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/Ped.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

class CalCalibSvc;

/** @class PedMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal pedestal calibration data
*/

class PedMgr : public CalibItemMgr {
 public:
  PedMgr() : 
    CalibItemMgr(CalibData::CAL_Ped),
    m_idealPeds(RngNum::N_VALS),
    m_idealPedSig(RngNum::N_VALS),
    m_idealCos(RngNum::N_VALS)
    {};

  /// get pedestal vals for given xtal/face/rng
  StatusCode getPed(CalXtalId xtalId,
                    float &avr,
                    float &sig,
                    float &cos);
 private:
  bool checkXtalId(CalXtalId xtalId) {
    if (!xtalId.validRange() || !xtalId.validFace())
      throw invalid_argument("Ped calib_type requires valid range & face info in CalXtalId."
                             " Programmer error");
    return true;
  }

  StatusCode loadIdealVals();

  LATWideIndex genIdx(CalXtalId xtalId) {return RngIdx(xtalId);}

  /// ped vals to use when calib db is down
  CalVec<RngNum,float> m_idealPeds;   
  /// ped sigma vals to use when calib db is down
  CalVec<RngNum,float> m_idealPedSig; 
  /// correlated ped cosine vals to use when calib db is down
  CalVec<RngNum,float> m_idealCos;    

  bool validateRangeBase(CalibData::RangeBase *rangeBase) {return true;}

};

#endif
