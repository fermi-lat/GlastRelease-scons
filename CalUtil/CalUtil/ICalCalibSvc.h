#ifndef ICalCalibSvc_H
#define ICalCalibSvc_H 1
// Include files
// Gaudi Generic
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Property.h"

// ROOT
#include "TSpline.h"

// Glast specific
#include "idents/CalXtalId.h"

// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_ICalCalibSvc("ICalCalibSvc", 0 , 1);

/*! @class ICalCalibSvc
 * \brief Abstract interface class for provision of GLAST LAT calorimiter calibration constants
 * \author Zach Fewtrell
 *
 * get*** functions are provided for each calibration type.  calibration constants are passed back by reference.
 *
 */

class ICalCalibSvc : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalCalibSvc; }

  /// retrieve CAL_ElecGain constant
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param gain output destination for gain constant
  /// \param sig output destination for sigma on gain
  virtual StatusCode getGain(const idents::CalXtalId &xtalId,
                             idents::CalXtalId::XtalFace face,
                             idents::CalXtalId::AdcRange range,
                             float &gain,
                             float &sig) = 0;


  /// retrieve integral non-linearity spline points
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param vals output const ref to vector of ADC values (y).
  /// \param dacs output const ref to vector of associated DAC values (x).
  /// \param error error values on ADC
  virtual StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const std::vector< float > *&vals,
                                  const std::vector< unsigned > *&dacs,
                                  float &error) = 0;

  /// retrieve integral non-linearity spline points as a spline object.
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param spline output const reference to ROOT TSpline3 object containing all int-nonlin points for given xtal/range.
  virtual StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const TSpline3 *&spline) = 0;


  /// retrieve light asymetry calibration constants
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param vals output vector of light asymetry constants, spacial resolution is one xtal width
  /// \param error output error value for vals
  virtual StatusCode getLightAsym(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const std::vector< float > *&vals,
                                  float &error) = 0;

  /// retrieve light attenuation calibration constant
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param att output light attenuation constant
  virtual StatusCode getLightAtt(const idents::CalXtalId &xtalId,
                                 idents::CalXtalId::XtalFace face,
                                 idents::CalXtalId::AdcRange range,
                                 float &att,
                                 float &norm) = 0;

  /// retrieve muon slope calibration constant
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param slope output muon slope constant
  /// \param error output sigma on slope
  virtual StatusCode getMuSlope(const idents::CalXtalId &xtalId,
                                idents::CalXtalId::XtalFace face,
                                idents::CalXtalId::AdcRange range,
                                float &slope,
                                float &error) = 0;

  /// retrieve pedestal calibration constant
  /// \param xtalId specify xtal log
  /// \param face specify xtal face
  /// \param range specify xtal range
  /// \param avr output pedestal
  /// \param sig output sigma on avr
  virtual StatusCode getPed(const idents::CalXtalId &xtalId,
                            idents::CalXtalId::XtalFace face,
                            idents::CalXtalId::AdcRange range,
                            float &avr,
                            float &sig,
                            float &cosAngle) = 0;

};

#endif // ICalCalibSvc_H
