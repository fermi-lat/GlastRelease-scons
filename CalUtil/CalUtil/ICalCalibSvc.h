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

/** @class ICalCalibSvc
 * @brief Interface class for CalCalibSvc
 *
 * Provides CalibData::CAL* constants in EZ interface.
 *
 * Attaches to one particular flavor of TDS data, more than one instance is possible for use of different flavors simulatenously.
 *
 * Is responsible for caching semi-static data such as IntNonlin and LightAsym spline objects to save recreating these objects for each interaction.  
 * Is also responsble for flushing the cache when the data becomes invalid.
 *
 *
 * Author:  Z. Fewtrell
 *
 */

class ICalCalibSvc : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalCalibSvc; }

  virtual StatusCode getGain(const idents::CalXtalId &xtalId,
                             idents::CalXtalId::XtalFace face,
                             idents::CalXtalId::AdcRange range,
                             float &gain,
                             float &sig) = 0;


  virtual StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const std::vector< float > *&vals,
                                  const std::vector< unsigned > *&dacs,
                                  float &error) = 0;

  virtual StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const TSpline3 *&intNonlinSpline) = 0;


  virtual StatusCode getLightAsym(const idents::CalXtalId &xtalId,
                                  idents::CalXtalId::XtalFace face,
                                  idents::CalXtalId::AdcRange range,
                                  const std::vector< float > *&vals,
                                  float &error) = 0;

  virtual StatusCode getLightAtt(const idents::CalXtalId &xtalId,
                                 idents::CalXtalId::XtalFace face,
                                 idents::CalXtalId::AdcRange range,
                                 float &att,
                                 float &norm) = 0;

  virtual StatusCode getMuSlope(const idents::CalXtalId &xtalId,
                                idents::CalXtalId::XtalFace face,
                                idents::CalXtalId::AdcRange range,
                                float &slope,
                                float &error) = 0;

  virtual StatusCode getPed(const idents::CalXtalId &xtalId,
                            idents::CalXtalId::XtalFace face,
                            idents::CalXtalId::AdcRange range,
                            float &avr,
                            float &sig,
                            float &cosAngle) = 0;

};

#endif // ICalCalibSvc_H
