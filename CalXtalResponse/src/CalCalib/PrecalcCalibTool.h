#ifndef PrecalcCalibTool_h
#define PrecalcCalibTool_h

// LOCAL INCLUDES
#include "CalXtalResponse/ICalCalibSvc.h"
#include "IPrecalcCalibTool.h"

// GLAST INCLUDES
#include "CalUtil/CalArray.h"
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IIncidentListener.h"

// STD INCLUDES

using namespace CalibData;
using namespace CalUtil;

/// \brief save time w/ certain pre-calculated values.
/// \note have to keep in sync w/ CalibDataSvc.
/// \note modeled after CalibItemMgr class
///  \note IIncidentListener allows me to check current validity for cached pre-calculated
///  calibration values.

class PrecalcCalibTool : public AlgTool,
                         virtual public IIncidentListener,
                         virtual public IPrecalcCalibTool
{
 public:
  PrecalcCalibTool(const string& type,
                   const string& name,
                   const IInterface* parent);

  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}

  /// return pedestal sigma converted to CIDAC scale
  StatusCode getPedSigCIDAC(RngIdx rngIdx, float &cidac);

  /// return trigger threshold in CIDAC scale
  StatusCode getTrigCIDAC(DiodeIdx diodeIdx, float &cidac);

  /// return trigger threshold in faceSignal (MeV) scale
  StatusCode getTrigMeV(DiodeIdx diodeIdx, float &mev);

  /// \brief return trigger threshold in proper adc range w/ associated range 
  /// & adc value
  ///
  /// this is needed bc trigger thresholds are often measured past the 
  /// saturation point of x8 range
  StatusCode getTrigRngADC(DiodeIdx diodeIdx, RngNum &rng, float &adc);
  
  /// return lac threshold in CIDAC scale
  StatusCode getLacCIDAC(FaceIdx faceIdx, float &lacCIDAC);


 private:
  /** \brief check calib validity period, (re)build local store if necessary
          
  needs to be called once per event (_before_ processing calibration data ;).
  Subsequent calls in same event will return immediately.
      
  */
  StatusCode updateCalib();          

  /// populate local data.
  StatusCode genLocalStore();

  /// serial # for current calibration source
  int m_serNoPed;      
  /// serial # for current calibration source
  int m_serNoINL;
  /// serial # for current calibration source
  int m_serNoTholdCI;    
  /// if false, current data has never been populated.
  bool m_initialized;
      
  bool m_isValid;

  /// trigger threholds in CIDAC scale
  CalArray<DiodeIdx, float> m_trigCIDAC;

  /// \brief trigger threholds in ADC scale (converted to proper x1 or x8 range)
  ///
  /// this is needed bc trigger thresholds are often measured past the 
  /// saturation point of x8 range
  CalArray<DiodeIdx, float> m_trigADC;

  /// \brief adc range for trigger threshold adc value
  CalArray<DiodeIdx, RngNum> m_trigRng;

  /// trigger threhold in faceSignal MeV
  CalArray<DiodeIdx, float> m_trigMeV;

  /// lac thresholds converted to CIDAC scale.
  CalArray<FaceIdx, float> m_lacCIDAC;
  /// pedestal noise sigma converted to CIDAC scale
  CalArray<RngIdx, float> m_pedSigCIDAC;

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      

  /// used to retrieve calib vals & check validity
  ICalCalibSvc *m_calCalibSvc;

  /// hook the BeginEvent so that we can check validity of XtalDigiPrecalc.
  void handle ( const Incident& inc ) {
    if ((inc.type() == "BeginEvent")) 
      m_isValid = false;
  }

  /// \brief Convert LEX8 adc value to LEX1 scale for given xtal face
  /// \note support LEX8 values > 4095 via extrapolation
  /// \param uldTholdX8 is needed for the calculation.  i know you already
  /// have it, so there's no sense in me retrieving it again.
  StatusCode lex8_to_lex1(FaceIdx faceIdx, 
                          float l8adc, float &l1adc);

  /// \brief Convert HEX8 adc value to HEX1 scale for given xtal face
  /// \note support HEX8 values > 4095 via extrapolation
  /// \param uldTholdX8 is needed for the calculation.  i know you already
  /// have it, so there's no sense in me retrieving it again.
  StatusCode hex8_to_hex1(FaceIdx faceIdx, 
                          float h8adc, float &h1adc);
};

#endif
