#ifndef ICalCalibSvc_H
#define ICalCalibSvc_H 1

// LOCAL INCLUDES

// GLAST INCLUDES
#include "idents/CalXtalId.h"
#include "CalibData/RangeBase.h" 

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Property.h"

// STD INCLUDES

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

using namespace std;
using namespace idents;

class ICalCalibSvc : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalCalibSvc; }

  /// retrieve MeVPerDac ratios for given xtal
  /// \param CalXtalId specify xtal log
  /// \param mpdLrg output MeVPerDac for large diode on both faces
  /// \param mpdSm output MeVPerDac for small diode on both faces
  virtual StatusCode getMeVPerDac(const CalXtalId &xtalId,
                                  CalibData::ValSig &mpdLrg,
                                  CalibData::ValSig &mpdSm) = 0;

  /// retrieve integral non-linearity vals for given xtal/face/rng
  /// \param CalXtalId specify xtal log
  /// \param vals output const ref to vector of ADC vals (y).
  /// \param dacs output const ref to vector of associated DAC vals (x).
  /// \param error error vals on ADC
  virtual StatusCode getIntNonlin(const CalXtalId &xtalId,
                                  const vector< float > *&adcs,
                                  const vector< unsigned > *&dacs,
                                  float &error) = 0;

  /// retrieve pedestal vals for given xtal/face/rng
  /// \param CalXtalId specify xtal log
  /// \param avr output ped (adc units)
  /// \param sig output sigma on avr
  virtual StatusCode getPed(const CalXtalId &xtalId,
                            float &avr,
                            float &sig,
                            float &cos) = 0;

  /// retrieve Asymmetry calibration information for one xtal
  /// \param CalXtalId specify xtal log
  /// \param large asymmetry array over position for large diode on both faces
  /// \param small asymmetry array over position for small diode on both faces
  /// \param MsmallPlrg asymmetry array over position for small diode on Minus face and lrg diode on Plus face
  /// \param PlrgMsmall asymmetry array over position for small diode on Plus face and lrg diode on Minus face
  /// \param xVals corresponding xvals for all of the accompanying asymmetry arrays
  virtual StatusCode getAsym(const CalXtalId &xtalId,
                             const vector<CalibData::ValSig> *&asymLrg,
                             const vector<CalibData::ValSig> *&asymSm,
                             const vector<CalibData::ValSig> *&asymNSPB,
                             const vector<CalibData::ValSig> *&asymPSNB,
                             const vector<float>  *&xVals) = 0;

  /// retrieve threshold calibration constants as measured w/ charnge injection
  /// \param CalXtalId specify xtal log
  /// \param FLE Fast Low Energy shaper threshold (Dac units)
  /// \param FHE Fast High Energy shaper threshold (Dac units)
  /// \param LAC Log accept threshold (Adc units)
  virtual StatusCode getTholdCI(const CalXtalId &xtalId,
                                CalibData::ValSig &FLE,
                                CalibData::ValSig &FHE,
                                CalibData::ValSig &LAC
                                ) = 0;

  /// retrieve Upper Level Discriminator threshold as measured w/ charnge injection for given xtal/face/rng
  /// \param CalXtalId specify xtal log
  /// \param output ULD threshold.
  virtual StatusCode getULDCI(const CalXtalId &xtalId,
                              CalibData::ValSig &ULDThold) = 0;

  /// retrieve pedestal calibration constants as measured during charge injection threshold testing.
  /// \param CalXtalId specify xtal log
  /// \param ped output ped val (Adc units)
  virtual StatusCode getPedCI(const CalXtalId &xtalId,
                              CalibData::ValSig &ped) = 0;

  /// retrieve threshold calibration constants as measured w/ muon calibration
  /// \param CalXtalId specify xtal log
  /// \param FLE Fast Low Energy shaper threshold (Dac units)
  /// \param FHE Fast High Energy shaper threshold (Dac units)
  virtual StatusCode getTholdMuon(const CalXtalId &xtalId,
                                  CalibData::ValSig &FLE,
                                  CalibData::ValSig &FHE
                                  ) = 0;

  /// retrieve pedestal calibration constants as measured during muon calibration threshold testing.
  /// \param CalXtalId specify xtal log
  /// \param ped output ped val (Adc units)
  virtual StatusCode getPedMuon(const CalXtalId &xtalId,
                                CalibData::ValSig &ped) = 0;

  virtual StatusCode evalDAC      (const CalXtalId &xtalId, double adc,   double &dac)  = 0;
  virtual StatusCode evalADC      (const CalXtalId &xtalId, double dac,   double &adc)  = 0;
  virtual StatusCode evalAsymLrg  (const CalXtalId &xtalId, double pos,   double &asym) = 0;
  virtual StatusCode evalPosLrg   (const CalXtalId &xtalId, double asym,  double &pos)  = 0;
  virtual StatusCode evalAsymSm   (const CalXtalId &xtalId, double pos,   double &asym) = 0;
  virtual StatusCode evalPosSm    (const CalXtalId &xtalId, double asym,  double &pos)  = 0;
  virtual StatusCode evalAsymNSPB (const CalXtalId &xtalId, double pos,   double &asym) = 0;
  virtual StatusCode evalPosNSPB  (const CalXtalId &xtalId, double asym,  double &pos)  = 0;
  virtual StatusCode evalAsymPSNB (const CalXtalId &xtalId, double pos,   double &asym) = 0;
  virtual StatusCode evalPosPSNB  (const CalXtalId &xtalId, double asym,  double &pos)  = 0;

};

#endif // ICalCalibSvc_H
