#ifndef IntNonlinMgr_H
#define IntNonlinMgr_H

// LOCAL
#include "CalibItemMgr.h"

// GLAST
#include "CalibData/Cal/IntNonlin.h"
#include "CalUtil/CalDefs.h"

// EXTLIB

// STD

using namespace CalDefs;
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
  StatusCode getIntNonlin(CalXtalId xtalId,
                          const vector< float > *&adcs,
                          const vector< float > *&dacs,
                          float &error);

  /// convert adc -> dac for given channel
  StatusCode evalDAC(CalXtalId xtalId, double adc, double &dac) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INL_SPLINE, xtalId, adc, dac);
  }
  
  /// convert dac -> adc for given channel
  StatusCode evalADC(CalXtalId xtalId, double dac, double &adc) {
    if (!checkXtalId(xtalId)) return StatusCode::FAILURE;
    return evalSpline(INV_INL_SPLINE, xtalId, dac, adc);
  }

 private:
  
  /// genearte locally stored spline functions & other data
  StatusCode genLocalStore();

  /// generate appropriate index for this calib_type
  LATWideIndex genIdx(CalXtalId xtalId) {return RngIdx(xtalId);}

  /// check that CalXtalId has all required fields populated
  bool checkXtalId(CalXtalId xtalId) {
    if (!xtalId.validRange() || !xtalId.validFace())
      throw invalid_argument("IntNonlin calib_type requires valid range "
                             "& face info in CalXtalId");
    return true;
  }

  /// enumerate spline function types for this calib_type
  enum SPLINE_TYPE {
    INL_SPLINE,
    INV_INL_SPLINE,
    N_SPLINE_TYPES
  };

  /// load ideal calibration values (in lieu of measured calibrations)
  StatusCode loadIdealVals();

  /// ideal ADC values for each energy range (shared by each channel)
  CalVec<RngNum, vector<float> >    m_idealADCs;
  /// ideal DAC values for each energy range (shared by each channel)
  CalVec<RngNum, vector<float> >    m_idealDACs;
  /// error value used throughout ideal intNonlin calibrations
  float                             m_idealErr;
  
  //-- LOCAL DATA STORE --//
  /** \brief Local copy of DAC values for each channel in consistent format

  TDS data can either store DAC values as integers or floats, either per channel or global, 
  depending on the version of the calibration file....  This allows me to keep a homogenous array
  of all floats, one set for every channel, regardless of the original data format.
   */
  CalVec<RngIdx, vector<float> > m_DACs;

  bool validateRangeBase(CalibData::RangeBase *rangeBase);
  
};

#endif
