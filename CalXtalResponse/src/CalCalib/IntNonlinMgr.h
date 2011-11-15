#ifndef IntNonlinMgr_H
#define IntNonlinMgr_H
// $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"
#include "CalibData/Cal/IntNonlin.h"

// EXTLIB

// STD
#include <memory>

class CalCalibSvc;

/** @class IntNonlinMgr
    @author Z.Fewtrell
    
    \brief Manage GLAST Cal integral non-linearity (CIDAC->ADC) calibration data.
*/

class IntNonlinMgr : public CalibItemMgr {
 public:
  IntNonlinMgr(CalCalibShared &ccsShared);

  /// get adc points for given cidac2adc curve
  const std::vector<float> *getInlAdc(const CalUtil::RngIdx rngIdx);

  /// get cidac points for given ciadc2adc curve
  const std::vector<float> *getInlCIDAC(const CalUtil::RngIdx rngIdx);

  /// evaluate cidac signal for given cidac2adc curve and adc value
  StatusCode evalCIDAC(const CalUtil::RngIdx rngIdx, const float adc, float &cidac);

  /// evaluate adc level for given cidac2adc curve and cidac value
  StatusCode evalADC(const CalUtil::RngIdx rngIdx, const float cidac, float &adc);

 private:
  
  /// generate locally stored spline functions & other data
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
  CalUtil::CalVec<CalUtil::RngIdx, std::vector<float> > m_CIDACs;

  /// ideal intNonlin calibration data (same for all xtals)
  CalUtil::CalArray<CalUtil::RngNum, std::auto_ptr<CalibData::IntNonlin> > m_idealINL;

  /// check calibration data for signle crystal
  bool validateRangeBase(CalibData::IntNonlin const*const intNonlin);
};

#endif
