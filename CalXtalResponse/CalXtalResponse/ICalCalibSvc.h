#ifndef ICalCalibSvc_H
#define ICalCalibSvc_H
//  $Header$

// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalibData/RangeBase.h"  // for ValSig
#include "CalUtil/CalArray.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"

// STD INCLUDES
#include <vector>

// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_ICalCalibSvc("ICalCalibSvc", 1, 0);

/*! @class ICalCalibSvc
 * \brief Abstract interface for provision of GLAST LAT calorimeter calib consts
 * \author Zach Fewtrell
 *
 * \note functions are provided for each calibration type.  
 calib consts are passed back by reference.
 *
 */

using namespace std;
using namespace idents;
using namespace CalUtil;

class ICalCalibSvc : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalCalibSvc; }


  /** \brief get MeVPerDac ratios for given xtal

  \param xtalIdx specify xtal log
  \param mpdLrg output MeVPerDac for large diode on both faces
  \param mpdSm output MeVPerDac for small diode on both faces
  */
  virtual StatusCode getMPD(XtalIdx xtalIdx,
                            CalibData::ValSig &mpdLrg,
                            CalibData::ValSig &mpdSm) = 0;
  

  /** \brief quick get both all MeVperDACvalues for single xtal
      \param xtalIdx specify cal xtal
      \param mpd output for 1 mevperdac per face
  */
  virtual StatusCode getMPD(XtalIdx xtalIdx, 
                            CalArray<DiodeNum, float> &mpd) = 0;


  /** \brief get integral non-linearity vals for given xtal, face/rng

  \param rngIdx specify xtal log, face, range
  \param vals output const ref to vector of ADC vals (y).
  \param dacs output const ref to vector of associated DAC vals (x).
  \param error error vals on ADC
  */
  virtual StatusCode getIntNonlin(RngIdx rngIdx,
                                  const vector< float > *&adcs,
                                  const vector< float > *&dacs,
                                  float &error) = 0;

  /** \brief get pedestal vals for given xtal, face/rng

  \param rngIdx specify xtal log, face, range
  \param avr output ped (adc units)
  \param sig output sigma on avr
  */
  virtual StatusCode getPed(RngIdx rngIdx,
                            float &avr,
                            float &sig,
                            float &cos) = 0;

  /** \brief quick get all pedestal & pedsig values for single xtal
      \param xtalIdx specify cal xtal
      \param ped one ped val per xtal/rng
      \param sig one pedestal sig per xtal/rng
  */
  virtual StatusCode getPed(XtalIdx xtalIdx,
                            CalArray<XtalRng, float> &peds,
                            CalArray<XtalRng, float> &sigs) = 0;

  /** \brief get Asymmetry calibration information for one xtal

  \param xtalIdx specify xtal log
  \param asymLrg asym array over position for large diode on both faces
  \param asymSm asym array over position for small diode on both faces
  \param asymNSPB asym array over position for small diode on
  Minus face and lrg diode on Plus face
  \param asymPSNB asym array over position for small diode on Plus
  face and lrg diode on Minus face

  \param xVals corresponding xvals for all of the accompanying asym arrays
  */
  virtual StatusCode getAsym(XtalIdx xtalIdx,
                             const vector<CalibData::ValSig> *&asymLrg,
                             const vector<CalibData::ValSig> *&asymSm,
                             const vector<CalibData::ValSig> *&asymNSPB,
                             const vector<CalibData::ValSig> *&asymPSNB,
                             const vector<float>  *&xVals) = 0;

  /** \brief get threshold calibration constants as measured w/ charge injection

  \param faceId specify xtal log, face
  \param FLE Fast Low Energy shaper threshold (Dac units)
  \param FHE Fast High Energy shaper threshold (Dac units)
  \param LAC Log accept threshold (Adc units)
  */
  virtual StatusCode getTholdCI(FaceIdx faceIdx,
                                CalibData::ValSig &FLE,
                                CalibData::ValSig &FHE,
                                CalibData::ValSig &LAC
                                ) = 0;

  /** \brief quick get all fle/fhe/lac vals for single xtal
      \param xtalIdx specify cal xtal
      \param trigThresh output trigger threshold for each xtal-diode
      \param lacThresh output LAC threshold for each xtal-face
  */

  virtual StatusCode getTholdCI(XtalIdx xtalIdx,
                                CalArray<XtalDiode, float> &trigThesh,
                                CalArray<FaceNum, float> &lacThresh) = 0;

  /** \brief get Upper Level Discriminator threshold as measured w/
      charge injection for given xtal, face/rng

      \param rngIdx specify xtal log, face, range
      \param output ULD threshold.
  */
  virtual StatusCode getULDCI(RngIdx rngIdx,
                              CalibData::ValSig &ULDThold) = 0;

  /** \brief quick get all ULD thesholds for single cal xtal
      \param
  */
  virtual StatusCode getULDCI(XtalIdx xtalIdx,
                              CalArray<XtalRng, float> &uldThold) = 0;

  /** \brief get pedestal calib constants from charge injection threshold tests.

  \param rngIdx specify xtal log, face, range
  \param ped output ped val (Adc units)
  */
  virtual StatusCode getPedCI(RngIdx rngIdx,
                              CalibData::ValSig &ped) = 0;

  /** \brief get threshold calibration constants as measured w/ muon calibration

  \param faceId specify xtal log, face
  \param FLE Fast Low Energy shaper threshold (Dac units)
  \param FHE Fast High Energy shaper threshold (Dac units)
  */
  virtual StatusCode getTholdMuon(FaceIdx faceIdx,
                                  CalibData::ValSig &FLE,
                                  CalibData::ValSig &FHE
                                  ) = 0;

  /** \brief get pedestal calib consts as measured during muon threshold testing.

  \param rngIdx specify xtal log, face, range
  \param ped output ped val (Adc units)
  */
  virtual StatusCode getPedMuon(RngIdx rngIdx,
                                CalibData::ValSig &ped) = 0;

  virtual StatusCode evalDAC      (RngIdx rngIdx, float adc,   float &dac)  
    = 0;
  virtual StatusCode evalADC      (RngIdx rngIdx, float dac,   float &adc)  
    = 0;
  virtual StatusCode evalAsymLrg  (XtalIdx xtalIdx, float pos,   float &asym) 
    = 0;
  virtual StatusCode evalPosLrg   (XtalIdx xtalIdx, float asym,  float &pos)  
    = 0;
  virtual StatusCode evalAsymSm   (XtalIdx xtalIdx, float pos,   float &asym) 
    = 0;
  virtual StatusCode evalPosSm    (XtalIdx xtalIdx, float asym,  float &pos)  
    = 0;
  virtual StatusCode evalAsymNSPB (XtalIdx xtalIdx, float pos,   float &asym) 
    = 0;
  virtual StatusCode evalPosNSPB  (XtalIdx xtalIdx, float asym,  float &pos)  
    = 0;
  virtual StatusCode evalAsymPSNB (XtalIdx xtalIdx, float pos,   float &asym) 
    = 0;
  virtual StatusCode evalPosPSNB  (XtalIdx xtalIdx, float asym,  float &pos)  
    = 0;

  /** \brief convert adc readout to MeV deposited at center of xtal

  \param rngIdx specify xtal log, face, range
  \param adcPed input ADC readout (pedestal subtracted)
  \param MeV output energy readout (in MeV
  */
  virtual StatusCode evalFaceSignal(RngIdx rngIdx, float adc, float &ene) 
    = 0;
};

#endif // ICalCalibSvc_H
