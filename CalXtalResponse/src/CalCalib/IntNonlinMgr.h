#ifndef IntNonlinMgr_H
#define IntNonlinMgr_H
// $Header$

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/IntNonlin.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB

// STD
#include <memory>

using namespace CalUtil;
using namespace idents;
using namespace CalibData;

class CalCalibSvc;

/** @class IntNonlinMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal integral non-linearity (CIDAC->ADC) calibration data.
*/

class IntNonlinMgr : public CalibItemMgr {
 public:
  IntNonlinMgr(CalCalibShared &ccsShared);

  const vector<float> *getInlAdc(CalUtil::RngIdx rngIdx);

  const vector<float> *getInlCIDAC(CalUtil::RngIdx rngIdx);

  StatusCode evalCIDAC(RngIdx rngIdx, float adc, float &cidac) {
    return evalSpline(INL_SPLINE, rngIdx, adc, cidac);
  }
  
  StatusCode evalADC(RngIdx rngIdx, float cidac, float &adc) {
    return evalSpline(INV_INL_SPLINE, rngIdx, cidac, adc);
  }

 private:
  
  /// genearte locally stored spline functions & other data
  StatusCode genLocalStore();

  /// enumerate spline function types for this calib_type
  enum SPLINE_TYPE {
    INL_SPLINE,
    INV_INL_SPLINE,
    N_SPLINE_TYPES
  };

  /// load ideal calibration values (in lieu of measured calibrations)
  StatusCode loadIdealVals();

  //-- LOCAL DATA STORE --//
  /** \brief Local copy of CIDAC values for each channel in consistent format

  TDS data can either store CIDAC values as integers or floats, either per channel or global, 
  depending on the version of the calibration file....  This allows me to keep a homogenous array
  of all floats, one set for every channel, regardless of the original data format.
  */
  CalArray<RngIdx, vector<float> > m_CIDACs;

  CalArray<RngNum, auto_ptr<IntNonlin> > m_idealINL;

  /// check ptr to TDS data
  bool validateRangeBase(IntNonlin *intNonlin);
};

#endif
