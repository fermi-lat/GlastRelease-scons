#ifndef IntNonlinMgr_H
#define IntNonlinMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/IntNonlin.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB

// STD

using namespace CalUtil;
using namespace idents;

class CalCalibSvc;

/** @class IntNonlinMgr
    @author Zachary Fewtrell
    
    \brief Manage GLAST Cal integral non-linearity (CIDAC->ADC) calibration data.
*/

class IntNonlinMgr : public CalibItemMgr {
 public:
  IntNonlinMgr();

  /// get integral non-linearity vals for given xtal/face/rng
  StatusCode getIntNonlin(RngIdx rngIdx,
                          const vector< float > *&adcs,
                          const vector< float > *&dacs,
                          float &error);

  /// convert adc -> dac for given channel
  StatusCode evalDAC(RngIdx rngIdx, float adc, float &dac) {
    return evalSpline(INL_SPLINE, rngIdx, adc, dac);
  }
  
  /// convert dac -> adc for given channel
  StatusCode evalADC(RngIdx rngIdx, float dac, float &adc) {
    return evalSpline(INV_INL_SPLINE, rngIdx, dac, adc);
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

  /// ideal ADC values for each energy range (shared by each channel)
  CalArray<RngNum, vector<float> >    m_idealADCs;
  /// ideal DAC values for each energy range (shared by each channel)
  CalArray<RngNum, vector<float> >    m_idealDACs;
  /// error value used throughout ideal intNonlin calibrations
  float                             m_idealErr;
  
  //-- LOCAL DATA STORE --//
  /** \brief Local copy of DAC values for each channel in consistent format

  TDS data can either store DAC values as integers or floats, either per channel or global, 
  depending on the version of the calibration file....  This allows me to keep a homogenous array
  of all floats, one set for every channel, regardless of the original data format.
  */
  CalArray<RngIdx, vector<float> > m_DACs;

  /// check ptr to TDS data
  bool validateRangeBase(CalibData::IntNonlin *intNonlin);
};

#endif
