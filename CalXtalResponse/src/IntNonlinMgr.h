#ifndef IntNonlinMgr_H
#define IntNonlinMgr_H 1

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/IntNonlin.h"
#include "CalUtil/CalDefs.h"

// EXTLIB

// STD

using namespace CalDefs;
using namespace idents;

class IntNonlinMgr : public CalibItemMgr {
 public:
  IntNonlinMgr(const IdealCalCalib &idealCalib);

  /// retrieve integral non-linearity vals for given xtal/face/rng
  StatusCode getIntNonlin(const CalXtalId &xtalId,
                          const vector< float > *&adcs,
                          const vector< unsigned > *&dacs,
                          float &error);
  StatusCode evalDAC(const CalXtalId &xtalId, double adc, double &dac) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INL_SPLINE, RngIdx(xtalId), adc, dac);
  }
  
  StatusCode evalADC(const CalXtalId &xtalId, double dac, double &adc) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_INL_SPLINE, RngIdx(xtalId), dac, adc);
  }

 private:
  CalibData::IntNonlin *getRangeBase(const CalXtalId &xtalId) {
    StatusCode sc = updateCache();
    if (sc.isFailure()) return NULL;
    return (CalibData::IntNonlin*)
      (m_rngBases[RngIdx(xtalId)]);
  }
  
  bool validateRangeBase(const CalXtalId &xtalId, 
                         CalibData::RangeBase *rngBase);
  
  StatusCode fillRangeBases();
  
  StatusCode genSplines();
  
  bool checkXtalId(const CalXtalId &xtalId) {
    if (!xtalId.validRange() || !xtalId.validFace())
      throw invalid_argument("IntNonlin calib_type requires valid range "
                             "& face info in CalXtalId");
    return true;
  }

  enum SPLINE_TYPE {
    INL_SPLINE,
    INV_INL_SPLINE,
    N_SPLINE_TYPES
  };

  StatusCode loadIdealVals();

  CalVec<RngNum, vector<float> >    m_idealADCs;
  CalVec<RngNum, vector<unsigned> > m_idealDACs;
  float                             m_idealErr;
  
};

#endif
