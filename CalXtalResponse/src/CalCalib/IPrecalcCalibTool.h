#ifndef IPrecalcCalibTool_h
#define IPrecalcCalibTool_h
// $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"


static const InterfaceID IID_IPrecalcCalibTool("IPrecalcCalibTool", 1, 0);


/** \brief Interface for tool which stores precacalculated calibration constants for GLAST Cal
*/
class IPrecalcCalibTool : virtual public IAlgTool {
 public:

  static const InterfaceID& interfaceID() { return IID_IPrecalcCalibTool; }

  virtual ~IPrecalcCalibTool() {}

  /// return pedestal sigma converted to CIDAC scale
  virtual StatusCode getPedSigCIDAC(CalUtil::RngIdx rngIdx, float &pedSigCIDAC) = 0;

  /// return trigger threshold in CIDAC scale
  virtual StatusCode getTrigCIDAC(CalUtil::DiodeIdx diodeIdx, float &trigCIDAC) = 0;
  
  /// return lac threshold in CIDAC scale
  virtual StatusCode getLacCIDAC(CalUtil::FaceIdx faceIdx, float &lacCIDAC) = 0;

  /// return trigger threshold in faceSignal (MeV) scale
  virtual StatusCode getTrigMeV(CalUtil::DiodeIdx diodeIdx, float &mev) = 0;

  /// \brief return trigger threshold in proper adc range w/ associated range 
  /// & adc value
  ///
  /// this is needed bc trigger thresholds are often measured past the 
  /// saturation point of x8 range
  virtual StatusCode getTrigRngADC(CalUtil::DiodeIdx diodeIdx, CalUtil::RngNum &rng, float &adc) = 0;


};

#endif
