#ifndef PedMgr_H
#define PedMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/Ped.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
// STD

using namespace CalDefs;
using namespace idents;

class PedMgr : public CalibItemMgr {
 public:
  PedMgr(const IdealCalCalib &idealCalib) : 
    CalibItemMgr(CalibData::CAL_Ped,
                 idealCalib),
    m_idealPeds(RngNum::N_VALS),
    m_idealPedSig(RngNum::N_VALS),
    m_idealCos(RngNum::N_VALS)
    {};

  /// retrieve pedestal vals for given xtal/face/rng
  StatusCode getPed(const CalXtalId &xtalId,
                    float &avr,
                    float &sig,
                    float &cos);
 private:
  CalibData::Ped *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::Ped*)(m_rngBases[RngIdx(xtalId)]);
  }
  
  bool validateRangeBase(const CalXtalId&, CalibData::RangeBase*)
    {return true;}
  
  StatusCode fillRangeBases();
  
  bool checkXtalId(const CalXtalId &xtalId) {
    if (!xtalId.validRange() || !xtalId.validFace())
      throw invalid_argument("Ped calib_type requires valid range & face info in CalXtalId");
    return true;
  }

  StatusCode loadIdealVals();

  CalVec<RngNum,float> m_idealPeds;   ///< ped vals to use when calib db is down
  CalVec<RngNum,float> m_idealPedSig; ///< ped sigma vals to use when calib db is down
  CalVec<RngNum,float> m_idealCos;    ///< correlated ped cosine vals to use when calib db is down
};

#endif
