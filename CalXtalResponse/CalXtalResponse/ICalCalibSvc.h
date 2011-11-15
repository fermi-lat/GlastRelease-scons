#ifndef ICalCalibSvc_H
#define ICalCalibSvc_H
//  $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES


// GLAST INCLUDES
#include "CalibData/RangeBase.h"  // for ValSig
#include "CalibData/Cal/Ped.h"
#include "CalibData/Cal/IntNonlin.h"
#include "CalibData/Cal/CalAsym.h"
#include "CalibData/Cal/CalMevPerDac.h"
#include "CalibData/Cal/CalTholdCI.h"
#include "CalibData/Cal/Xpos.h"
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"

// STD INCLUDES

// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_ICalCalibSvc("ICalCalibSvc", 2, 2);

/*! @class ICalCalibSvc
 * \brief Abstract interface for provision of GLAST LAT calorimeter calib constants
 * \author Z.Fewtrell
 *
 *
 */
class ICalCalibSvc : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalCalibSvc; }

  virtual ~ICalCalibSvc() {};


  /** \brief get MeVPerDac ratios for given xtal

  \param xtalIdx specify xtal log
  \return null on error. 
  */
  virtual const CalibData::CalMevPerDac *getMPD(CalUtil::XtalIdx xtalIdx) = 0;


  /** \brief get adc value list for intNonlin (adc2cidac) curves.

      \note since internal storage of adc & cidac information has changed, we are providing
      this generic function in leiu of the actual internal IntNonlin object.

      \return 0 on error
   */
  virtual const std::vector<float> *getInlAdc(CalUtil::RngIdx rngIdx) = 0;

  /** \brief get cidac value list for intNonlin (adc2cidac) curves.

      \note since internal storage of adc & cidac information has changed, we are providing
      this generic function in leiu of the actual internal IntNonlin object.

      \return 0 on error
   */
  virtual const std::vector<float> *getInlCIDAC(CalUtil::RngIdx rngIdx) = 0;

  /** \brief get pedestal vals for given xtal, face/rng

  \param rngIdx specify xtal log, face, range
  \return null on error.
  */
  //  virtual const CalibData::Ped *getPed(CalUtil::RngIdx rngIdx) = 0;
  virtual StatusCode getPed(CalUtil::RngIdx rngIdx, float &ped)=0;
  virtual StatusCode getPedSig(CalUtil::RngIdx rngIdx, float &sig)=0;


  /** \brief get Asymmetry calibration information for one xtal

  \param xtalIdx specify xtal log
  \return 0 on error.

  */
  virtual const CalibData::CalAsym *getAsym(CalUtil::XtalIdx xtalIdx) = 0;
  
  /** \brief get global xtal position list assocated w/ Asymmetry 
      arrays in CalAsym objects
      \return 0 on error
  */
  virtual const CalibData::Xpos *getAsymXpos() = 0;


  /** \brief get threshold calibration constants as measured w/ charge injection

  \param faceId specify xtal log, face
  \return 0 on error
  */
  virtual const CalibData::CalTholdCI *getTholdCI(CalUtil::FaceIdx faceIdx) = 0;

  /// convert adc -> cicidac for given channel
  virtual StatusCode evalCIDAC (CalUtil::RngIdx rngIdx, float adc,   float &cidac)  
    = 0;

  /// convert cicidac -> adc for given channel
  virtual StatusCode evalADC (CalUtil::RngIdx rngIdx, float cidac,   float &adc)  
    = 0;

  /// \brief convert xtal longitudinal pos -> signal asymmetry for given xtal & 
  /// diode pair
  virtual StatusCode evalAsym(CalUtil::XtalIdx xtalIdx, 
                              CalUtil::AsymType asymType, 
                              float pos,   float &asym) = 0;

  /// convert xtal signal asymmetry -> longitudinal pos for given xtal & diodes
  virtual StatusCode evalPos (CalUtil::XtalIdx xtalIdx, 
                              CalUtil::AsymType asymType, 
                              float asym,  float &pos)  = 0;

  /// evaluate optical signal asymmetry for deposit @ center of xtal
  virtual StatusCode getAsymCtr(CalUtil::XtalIdx xtalIdx, 
                                CalUtil::AsymType asymType, 
                                float &asymCtr) = 0;                                

  /** \brief convert adc readout to MeV deposited at center of xtal

  \param rngIdx specify xtal log, face, range
  \param adcPed input ADC readout (pedestal subtracted)
  \param MeV output energy readout (in MeV
  */
  virtual StatusCode evalFaceSignal(CalUtil::RngIdx rngIdx, 
                                    float adc, 
                                    float &ene) = 0;

  /** \brief get ratio of MeV/CIDAC for given diode where MeV is energy deposited @ center of xtal
      \note this calculation should include comination of overall mevPerDAC & asymmetry calibrations.
  */
  virtual StatusCode getMPDDiode(CalUtil::DiodeIdx diodeIdx, float &mpdDiode) = 0;

  /// retrieve serial # for current pedestal calibration data
  virtual int getSerNoPed() = 0;
  /// retrieve serial # for current pedestal intNonlin data
  virtual int getSerNoINL() = 0;
  /// retrieve serial # for current pedestal asymmetry data
  virtual int getSerNoAsym() = 0;
  /// retrieve serial # for current pedestal mevPerDAC data
  virtual int getSerNoMPD() = 0;
  /// retrieve serial # for current pedestal tholdCI data
  virtual int getSerNoTholdCI() = 0;

};

#endif // ICalCalibSvc_H
